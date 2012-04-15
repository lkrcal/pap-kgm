/**
 *  \file kgm.cpp
 *  \brief KGM
 *  \author krcallub, varkoale
 *  \version 0.1
 *
 *  --
 */

#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/timer.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/program_options.hpp>

#include <algorithm>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <utility>
#include <fstream>
#include <stack>
#include <queue>

#include <stdint.h>
#include <omp.h>

using namespace std;
using namespace boost;

struct kgm_vertex_properties {
    kgm_vertex_properties() :
        state(false),
        degree(0) {}
    bool state;
    uint16_t degree;
};

enum VERTEX_KGM_STATE {
    NEW,
    OPEN,
    __VERTEX_KGM_STATE_SIZE
};

typedef adjacency_list<
                listS, // std::list
                vecS, // std::vector
                undirectedS,
                kgm_vertex_properties> ugraph;

typedef graph_traits<ugraph>::edge_descriptor       kgm_edge_descriptor;
typedef graph_traits<ugraph>::vertex_descriptor     kgm_vertex_descriptor;
typedef graph_traits<ugraph>::edge_iterator         kgm_edge_iterator;
typedef graph_traits<ugraph>::vertex_iterator       kgm_vertex_iterator;
typedef graph_traits<ugraph>::adjacency_iterator    kgm_adjacency_iterator;

struct dfs_state {
    kgm_vertex_iterator v_it;
    kgm_vertex_iterator v_it_end;
    kgm_adjacency_iterator a_it;
    kgm_adjacency_iterator a_it_end;
};
struct i_dfs_state {
    uint16_t v;
    uint16_t a;
};
typedef std::vector<dfs_state> dfs_stack;
typedef std::vector<i_dfs_state> i_dfs_stack;
typedef std::vector<uint16_t> degree_stack;

struct kgm_task_struct
{
	std::vector<i_dfs_state> inactiveDfsStack;
	uint16_t newVItStart;
	uint16_t newVItEnd;
	uint16_t degree;
};

typedef boost::shared_ptr<kgm_task_struct> kgm_task;
typedef std::queue<kgm_task> kgm_task_queue;

ostream& operator<< (ostream& out, const kgm_adjacency_iterator& ai);
ostream& operator<< (ostream& out, const kgm_vertex_iterator& ai);
ostream& operator<< (ostream& out, const dfs_state& state);
ostream& operator<< (ostream& out, const std::pair<dfs_state,ugraph>& state);
ostream& operator<< (ostream& out, const dfs_stack& stack);


// RUNTIME -----------------------------------------------------

//ugraph g; // TODO private ( TODO needs to be copied via boost::graph::copy.hpp,
		  // cannot be made private by openmp's #pragma private)
uint32_t KGM_GRAPH_SIZE; // shared, RO
uint16_t KGM_UPPER_BOUND = 30; // shared // TODO add atomic to Write manimulation

const uint32_t KGM_LOWER_BOUND = 2;
const uint32_t KGM_START_NODE = 0;
const uint64_t KGM_REPORT_INTERVAL = 0x10000000;
const uint32_t KGM_GIVEAWAY_INTERVAL = 0x100;

boost::scoped_ptr<boost::timer> KGM_TIMER; // shared

bool running = true; // TODO private
uint16_t stoppedThreads = 0; // TODO shared

kgm_task_queue KGM_TASK_QUEUE;

// end (RUNTIME) -----------------------------------------------


int openMpThreadNumber()
{
	return omp_get_thread_num();
}
int openMpTotalThreads()
{
	return omp_get_num_threads();
}
// TODO Optimal solution found
void finish();

bool isValid_dfs_state(const dfs_state& dfsState, const ugraph& graph)
{
    if (dfsState.v_it == dfsState.v_it_end)
        return false;
    if (dfsState.a_it == dfsState.a_it_end)
        return false;
    if (graph[*(dfsState.v_it)].state == false ||
        graph[*(dfsState.a_it)].state == true)
        return false;
    return true;
}

bool iterate_dfs_state(dfs_state& dfsState, const ugraph& graph)
{
    if (graph[*(dfsState.v_it)].state == true)
    {
        for (++dfsState.a_it; dfsState.a_it != dfsState.a_it_end; ++dfsState.a_it)
        {
            if (graph[*(dfsState.a_it)].state == false)
                return true;
        }
    }
    for (++dfsState.v_it; dfsState.v_it != dfsState.v_it_end; ++dfsState.v_it)
    {
        if (graph[*(dfsState.v_it)].state == false)
            continue;

        for(tie(dfsState.a_it, dfsState.a_it_end) = adjacent_vertices(*(dfsState.v_it),graph);
                dfsState.a_it != dfsState.a_it_end; ++dfsState.a_it)
        {
            if (graph[*(dfsState.a_it)].state == false)
                return true;
        }
    }
    return false;
}

bool create_dfs_state(dfs_state& newState, const ugraph& graph)
{
    for (tie(newState.v_it, newState.v_it_end) = vertices(graph);
            newState.v_it != newState.v_it_end; ++newState.v_it) // nodes that are part of the skeleton
    {
        if (graph[*(newState.v_it)].state == false)
            continue;

        for(tie(newState.a_it, newState.a_it_end) = adjacent_vertices(*(newState.v_it),graph); // neighbours
                newState.a_it != newState.a_it_end; ++newState.a_it)
        {
            if (graph[*(newState.a_it)].state == false)
                return true;
        }
    }
    return false;
}

ostream& operator<< (ostream& out, const kgm_adjacency_iterator& ai)
{
    out << *ai;
    return out;
}

ostream& operator<< (ostream& out, const kgm_vertex_iterator& ai)
{
    out << *ai;
    return out;
}

ostream& operator<< (ostream& out, const dfs_state& state)
{
    out << "[" << *(state.v_it) << "->" << *(state.a_it) << "]";
    return out;
}

ostream& operator<< (ostream& out, const std::pair<dfs_state,ugraph>& state)
{
    out << "[";
    dfs_state tmp (state.first);
    bool first = true;
    do {
        if (first)
            first = false;
        else
            out << ", ";
        out << *(tmp.v_it) << "->" << *(tmp.a_it);
    }
    while (iterate_dfs_state(tmp, state.second));
    out << "]";
    return out;
}

ostream& operator<< (ostream& out, const dfs_stack& stack)
{
    for (dfs_stack::const_iterator it = stack.begin();
            it != stack.end(); ++it)
        out << *it;
    return out;
}
ostream& operator<< (ostream& out, const i_dfs_stack& stack)
{
    for (i_dfs_stack::const_iterator it = stack.begin();
            it != stack.end(); ++it)
        out << "[" << (*it).v << "--" << (*it).a << "]";
    return out;
}

void dfs_step(
		i_dfs_stack& inactiveDfsStack,
        dfs_stack& stack,
        degree_stack& degStack,
        ugraph& graph,
        uint16_t& degreeLimit) // so far the lowest skeleton degree
{
    if (stack.empty())
        return;

    if (degreeLimit <= KGM_LOWER_BOUND){
        stack.clear();
        return;
    }

    if (graph[*(stack.back().a_it)].state == false)
    {
        graph[*(stack.back().v_it)].degree += 1;
        graph[*(stack.back().a_it)].degree = 1;
        graph[*(stack.back().a_it)].state = true;

        degStack.push_back(std::max(graph[*(stack.back().v_it)].degree, degStack.back()));

        if (degStack.back() >= degreeLimit)
        {   // not going to be a better solution
            degStack.pop_back();
            graph[*(stack.back().v_it)].degree -= 1;
            graph[*(stack.back().a_it)].degree = 0;
            graph[*(stack.back().a_it)].state = false;

            if (!iterate_dfs_state(stack.back(), graph))
                stack.pop_back();
            return;
        }

        dfs_state newState;
        if (create_dfs_state(newState, graph))
            stack.push_back(newState);
        else
        {
            degreeLimit = degStack.back();
            // TODO possible spanning tree improvement
            std::cout << openMpThreadNumber() << ": New spanning tree of degree: "<< degreeLimit << std::endl
                    << openMpThreadNumber() << ": " << inactiveDfsStack << stack << std::endl
                    << openMpThreadNumber() << ": of total size: " << inactiveDfsStack.size() + stack.size() << std::endl;


			#pragma omp critical(kgmUpperBound)
			{
            	if (degreeLimit < KGM_UPPER_BOUND)
            		KGM_UPPER_BOUND = degreeLimit;
			}
            if(degreeLimit == KGM_LOWER_BOUND)
            {
            	// Terminate all threads
				#pragma omp critical(stoppedThreads)
            	{
            		stoppedThreads = openMpTotalThreads();
            	}
            }
        }

        return;
    }

    dfs_state prev (stack.back());
    graph[*(prev.v_it)].degree -= 1;
    graph[*(prev.a_it)].degree = 0;
    graph[*(prev.a_it)].state = false;

    if (!degStack.empty())
    degStack.pop_back();

    if (!stack.empty())
        if (!iterate_dfs_state(stack.back(), graph)) // not able to iterate further, branch is removed from the stack
            stack.pop_back();
}

void readGraphFromFile(std::string filename, ugraph& g) {
	std::ifstream in(filename.c_str());
	std::stringstream buffer;
	buffer << in.rdbuf();
	std::string contents(buffer.str());

	std::vector<std::string> lines;
	boost::split(lines, contents, boost::is_any_of("\n"));

	int m = lines.size()-2, n = lines[lines.size()-2].length();
	std::cout << "loaded graph of size: " << m << "*" << n << std::endl;
	KGM_GRAPH_SIZE = n;
	if (m != n && m <= 1)
	{
		std::cerr << "Input data error" << std::endl;
		exit(-3);
	}

	g = ugraph(KGM_GRAPH_SIZE);
	std::vector<std::string>::iterator it = lines.begin();
	++it;
	for (uint32_t i = 0; it != lines.end(); ++it, ++i)
	{
		if ((*it).empty())
			continue;
		for (uint32_t j = 0; j <= i; ++j)
		{
			if ((*it)[j] == '1')
				add_edge(i,j,g);
		}
	}
}

bool hasExtraWork(dfs_stack& dfsStack) {
    if (dfsStack.size() <= 1)
        return false;
    for (dfs_stack::iterator it = dfsStack.begin(); it < (dfsStack.end()-1); ++it)
    {
        if ((*it).v_it != (*it).v_it_end)
            return true;
    }
	return false;
}

bool giveWork(ugraph& g, i_dfs_stack& inactiveDfsStack, dfs_stack& dfsStack, degree_stack& degreeStack)
{
	int difference, diffRatio;
	kgm_vertex_iterator newVIt, newVIt_end;
	dfs_stack::iterator it;
	bool ok = false;

	kgm_task newTask;
	newTask.reset(new kgm_task_struct);

	// Cuts off last state from active stack
	for (it = dfsStack.begin(); it < (dfsStack.end()-1); ++it)
	{
		difference = (*it).v_it_end-(*it).v_it;
		if (difference > 1)
		{
			diffRatio = ((*it).v_it_end-(*it).v_it+1)/2;
			newVIt = (*it).v_it+diffRatio;
			newVIt_end = (*it).v_it_end;
			ok = true;
			break;
		}
	}
	if (!ok)
		return false;
    (*it).v_it_end = newVIt;

    newTask->newVItStart = (uint16_t)(*newVIt);
    newTask->newVItEnd = (uint16_t)(*newVIt_end);

    newTask->degree = degreeStack.at(it-dfsStack.begin());

    dfs_stack::iterator git;
	i_dfs_stack::iterator igit;
	for (igit = inactiveDfsStack.begin(); igit != inactiveDfsStack.end(); ++igit)
	{
		i_dfs_state iDfsState;
		iDfsState.v = (*igit).v;
		iDfsState.a = (*igit).a;
		newTask->inactiveDfsStack.push_back(iDfsState);
	}
	for (git = dfsStack.begin(); git != it; ++git)
	{
		i_dfs_state iDfsState;
		iDfsState.v = (uint16_t)(*((*git).v_it));
		iDfsState.a = (uint16_t)(*((*git).a_it));
		newTask->inactiveDfsStack.push_back(iDfsState);
	}

	#pragma omp critical(kgmtasklist)
	{
		KGM_TASK_QUEUE.push(newTask);
	}
	return true;
}

bool grabWork(ugraph& g, i_dfs_stack& inactiveDfsStack, dfs_stack& dfsStack, degree_stack& degreeStack)
{
	kgm_task grabbedTask;
	bool ok;
	#pragma omp critical(kgmtasklist)
	{
		if (KGM_TASK_QUEUE.empty())
			ok = false;
		else
		{
			grabbedTask = KGM_TASK_QUEUE.front();
			KGM_TASK_QUEUE.pop();
		}
	}
	if (!ok)
		return false;

    // Clear degree stack, put first state's degree into the stack
	degreeStack.clear();
	degreeStack.push_back(grabbedTask->degree);

	// Clear graph state
	kgm_vertex_iterator vi, vi_end;
	tie(vi, vi_end) = vertices(g);
	g[(*vi)].state = true;
	g[(*vi)].degree = 0;
	++vi;
	for (; vi != vi_end; ++vi)
	{
		g[(*vi)].state = false;
		g[(*vi)].degree = 0;
	}
	inactiveDfsStack.clear();

	// Initialize graph state from task
	for (i_dfs_stack::iterator it = grabbedTask->inactiveDfsStack.begin();
				it != grabbedTask->inactiveDfsStack.end(); ++it)
	{
		i_dfs_state istate;
		istate.v = (*it).v;
		istate.a = (*it).a;
		inactiveDfsStack.push_back(istate);
		g[istate.v].state = true;
		g[istate.v].degree++;
		g[istate.a].state = true;
		g[istate.a].degree++;
	}

	dfsStack.clear();
	dfs_state firstState;
	tie(firstState.v_it, firstState.v_it_end) = vertices(g);
	firstState.v_it += grabbedTask->newVItStart;
	firstState.v_it_end = firstState.v_it + (grabbedTask->newVItEnd - grabbedTask->newVItStart);
	std::cout << openMpThreadNumber() << ": firstState.v_it_end = firstState.v_it + "
			<< grabbedTask->newVItEnd-grabbedTask->newVItStart << std::endl;
	tie(firstState.a_it, firstState.a_it_end) = adjacent_vertices(*(firstState.v_it),g);

	std::cout << openMpThreadNumber() << ": raw first state received: " << firstState << std::endl;

	if (!isValid_dfs_state(firstState, g))
	{
		std::cout << openMpThreadNumber() << ": is invalid" << std::endl;
		iterate_dfs_state(firstState, g);
		std::cout << openMpThreadNumber() << ": iterated to: " << firstState << std::endl;
	} else {
		std::cout << openMpThreadNumber() << ": is valid" << std::endl;
	}
	if (!isValid_dfs_state(firstState, g))
	{
		std::cout << openMpThreadNumber() << ": is invalid after iteration" << std::endl
				  << openMpThreadNumber() << ": has accepted no valid work!" << std::endl;
	}
	else
	{
		std::cout << openMpThreadNumber() << ": successfully accepted work" << std::endl
				<< "   dfs stack size: " << dfsStack.size() << std::endl
				<< "   degree: " << degreeStack.back() << std::endl
				<< "   edges before this state: " << inactiveDfsStack.size() << std::endl;
		dfsStack.push_back(firstState);
	}
	return true;
}


void divideWork(ugraph& g, i_dfs_stack& inactiveDfsStack, dfs_stack& dfsStack, degree_stack& degreeStack) {
    std::cout << ": divideWork()" << std::endl;
	int threadNum = 0;
	while(hasExtraWork(dfsStack) && threadNum < openMpTotalThreads()) {
		giveWork(g, inactiveDfsStack, dfsStack, degreeStack);
		++threadNum;
	}
}

void initStack(ugraph& g, i_dfs_stack& inactiveDfsStack, dfs_stack& dfsStack, degree_stack& degreeStack) {
	dfs_state firstState;
	g[KGM_START_NODE].state = true;
	if (!create_dfs_state(firstState,g))
	{
		std::cerr << "Failed to initialize first kgm_process_state" << std::endl;
		exit(-4);
	}
	dfsStack.push_back(firstState);
	degreeStack.push_back(0);
}


void iterateStack(ugraph& g, i_dfs_stack& inactiveDfsStack, dfs_stack& dfsStack, degree_stack& degreeStack) {

	uint64_t KGM_REPORT_NEXT = KGM_REPORT_INTERVAL;
	uint64_t KGM_STEPS = 0;

	while(running) {
		while (!dfsStack.empty())
		{
			++KGM_STEPS;

			if (dfsStack.empty())
				break;


			dfs_step(inactiveDfsStack, dfsStack, degreeStack, g, KGM_UPPER_BOUND);

			if (KGM_STEPS % KGM_GIVEAWAY_INTERVAL == 0)
			{
				std::cout << openMpThreadNumber() << ": " << "Attepmting to give away work..." << std::endl;
				if (giveWork(g, inactiveDfsStack, dfsStack, degreeStack))
				{
					std::cout << openMpThreadNumber() << ": " << "Attepmt to give away work successful" << std::endl;
				}
				else
				{
					std::cout << openMpThreadNumber() << ": " << "Attempt to give away work successful" << std::endl;
				}

			}

			if (KGM_STEPS >= KGM_REPORT_NEXT)
			{
				std::cout << ": Stack report at " << KGM_TIMER->elapsed() << ":" << std::endl
						<< inactiveDfsStack << dfsStack << std::endl;
				KGM_REPORT_NEXT += KGM_REPORT_INTERVAL;
			}
		}

		if (stoppedThreads >= openMpTotalThreads())
		{
			std::cout << openMpThreadNumber() << ": " << "finished (optimal solution found)" << std::endl;
			running = false;
			break;
		}

		// TODO end on first failure is too broadminded
		if (!grabWork(g, inactiveDfsStack, dfsStack, degreeStack))
		{
			std::cout << openMpThreadNumber() << ": " << "finished" << std::endl;
			running = false;
			#pragma omp atomic
			stoppedThreads++;
		}
	}

}

int main(int argc, char ** argv) {

    if (argc <= 1)
    {
        std::cerr << "Not enough arguments." << std::endl;
        return -1;
    }

    std::string filename (argv[1]);
	if (filename.empty())
		return -1;
	ugraph prototype;
	readGraphFromFile(filename, prototype);


	i_dfs_stack inactiveDfsStack;
	dfs_stack dfsStack;
	degree_stack degreeStack;
	ugraph g; // TODO copy
    initStack(prototype, inactiveDfsStack, dfsStack, degreeStack);

    // TODO run for some time (expand)

    // TODO divide work

    KGM_TIMER.reset(new boost::timer);

    iterateStack(prototype, inactiveDfsStack, dfsStack, degreeStack);

	std::cout << ": ***** ENDED in time " << KGM_TIMER->elapsed() << std::endl;

	return 0;
}

