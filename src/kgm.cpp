/**
 *  \file kgm.cpp
 *  \brief KGM
 *  \author krcallub, kopecji5, varkoale
 *  \version 0.1
 *
 *  --
 */

#define _POSIX_C_SOURCE 199309L

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

#include <stddef.h>
#include <stdint.h>
#include <omp.h>
#include <time.h>
#include <ctime>

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
	uint32_t id;
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

uint32_t KGM_GRAPH_SIZE;
uint16_t KGM_UPPER_BOUND = 30;

const uint32_t KGM_LOWER_BOUND = 2;
const uint32_t KGM_START_NODE = 0;
const uint64_t KGM_REPORT_INTERVAL = 0x10000000;
const uint32_t KGM_GIVEAWAY_INTERVAL = 0x10000;

uint32_t NUMBER_OF_THREADS;

boost::scoped_ptr<boost::timer> KGM_TIMER;

uint16_t stoppedThreads = 0;

uint32_t KGM_TASK_COUNTER = 0;
kgm_task_queue KGM_TASK_QUEUE;

uint16_t finalDegree = KGM_UPPER_BOUND;
i_dfs_stack finalInactiveDfsStack;
dfs_stack finalDfsStack;

// end (RUNTIME) -----------------------------------------------


inline int openMpThreadNumber(){return omp_get_thread_num();}

inline int openMpTotalThreads(){return omp_get_num_threads();}

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
        uint16_t& degreeLimit) // so far the lowest spanning tree degree
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
            std::cout << openMpThreadNumber() << ": New spanning tree of degree: "<< degreeLimit << std::endl
                    << openMpThreadNumber() << ": " << inactiveDfsStack << stack << std::endl
                    << openMpThreadNumber() << ": of total size: " << inactiveDfsStack.size() + stack.size() << std::endl;

			#pragma omp critical(newSpanningTree)
            {
            	finalDegree = degreeLimit;
            	finalInactiveDfsStack = inactiveDfsStack;
            	finalDfsStack = stack;
            }

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
	#pragma omp critical(kgmTaskCounter)
	{
		newTask->id = KGM_TASK_COUNTER++;
	}

	#ifdef VERBOSE
		std::cout << openMpThreadNumber() << ": Giving away task " << newTask->id << std::endl;
	#endif

	#pragma omp critical(kgmtasklist)
	{
		KGM_TASK_QUEUE.push(newTask);
	}
	return true;
}

bool grabWork(ugraph& g, i_dfs_stack& inactiveDfsStack, dfs_stack& dfsStack, degree_stack& degreeStack)
{
	kgm_task grabbedTask;
	bool ok = true;
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

	#ifdef VERBOSE
		std::cout << openMpThreadNumber() << ": Accepted task " << grabbedTask->id << std::endl;
	#endif
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
	tie(firstState.a_it, firstState.a_it_end) = adjacent_vertices(*(firstState.v_it),g);

	if (!isValid_dfs_state(firstState, g))
		iterate_dfs_state(firstState, g);
	if (isValid_dfs_state(firstState, g))
		dfsStack.push_back(firstState);
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


void iterateStack(ugraph& g, i_dfs_stack& inactiveDfsStack, dfs_stack& dfsStack, degree_stack& degreeStack, bool& running) {

	uint64_t KGM_REPORT_NEXT = KGM_REPORT_INTERVAL;
	uint64_t KGM_GIVEAWAY_NEXT = KGM_GIVEAWAY_INTERVAL;
	uint64_t KGM_STEPS = 0;

	while(true) {
		while (!dfsStack.empty())
		{
			++KGM_STEPS;

			dfs_step(inactiveDfsStack, dfsStack, degreeStack, g, KGM_UPPER_BOUND);

			if (KGM_STEPS >= KGM_GIVEAWAY_NEXT)
			{
				giveWork(g, inactiveDfsStack, dfsStack, degreeStack);
				KGM_GIVEAWAY_NEXT += KGM_GIVEAWAY_INTERVAL;
			}


			if (KGM_STEPS >= KGM_REPORT_NEXT)
			{
				std::cout << openMpThreadNumber() << ": Stack report at " << KGM_TIMER->elapsed() << ":" << std::endl
						<< inactiveDfsStack << dfsStack << std::endl;
				KGM_REPORT_NEXT += KGM_REPORT_INTERVAL;
			}
		}

		if (stoppedThreads >= openMpTotalThreads())
		{
			std::cout << openMpThreadNumber() << ": " << "Finished!" << std::endl;
			running = false;
			break;
		}

		if (!grabWork(g, inactiveDfsStack, dfsStack, degreeStack))
		{
			if (running)
			{
				running = false;
				#pragma omp atomic
					stoppedThreads++;
			}
			usleep(1000);
		}
		else
		{
			if (!running)
			{
				running = true;
				#pragma omp atomic
					stoppedThreads--;
			}
		}
	}

}

int main(int argc, char ** argv) {

    if (argc <= 2)
    {
        std::cerr << "Not enough arguments." << std::endl;
        return -1;
    }
    istringstream(argv[2]) >> NUMBER_OF_THREADS;
    if (NUMBER_OF_THREADS <= 0  || NUMBER_OF_THREADS > 256)
    {
        std::cerr << "Number of threads is out of <1,256>." << std::endl;
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

    KGM_TIMER.reset(new boost::timer);

    timespec startWTime, endWTime;
    long int totalWTime;
    int totalWTimeSec;


    clock_gettime(CLOCK_MONOTONIC, &startWTime);

	#pragma omp parallel \
		private(inactiveDfsStack,dfsStack,degreeStack) \
		shared(KGM_TASK_COUNTER, KGM_TASK_QUEUE,KGM_UPPER_BOUND,stoppedThreads) \
		num_threads(NUMBER_OF_THREADS)
    {

		ugraph g;
		#pragma omp critical(ugraphCopy)
    	{
			g = prototype;
    	}
    	bool running = true;

		#pragma omp master
    	{
    		initStack(g, inactiveDfsStack, dfsStack, degreeStack);
    	}
		iterateStack(g, inactiveDfsStack, dfsStack, degreeStack, running);
    }

    clock_gettime(CLOCK_MONOTONIC, &endWTime);
    totalWTime = (endWTime.tv_nsec-startWTime.tv_nsec);
    totalWTimeSec = endWTime.tv_sec-startWTime.tv_sec;

	std::cout << ": *****************************************" << std::endl;
	std::cout << ": ***** CPU TIME (s): " << KGM_TIMER->elapsed() << std::endl;
	std::cout << ": ***** REAL TIME (s): " << totalWTimeSec << std::endl;
	std::cout << ": ***** REAL TIME (ns): " << totalWTime << std::endl;
	std::cout << ": ***** THREADS: " << NUMBER_OF_THREADS << std::endl;
	std::cout << ": ***** DEGREE: " << finalDegree << std::endl;
	std::cout << ": ***** SPANNING TREE: " << finalInactiveDfsStack << finalDfsStack << std::endl;
	std::cout << ": *****************************************" << std::endl;
	return 0;
}

