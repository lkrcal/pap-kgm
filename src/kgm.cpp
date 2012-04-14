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

#include <stdint.h>

#include <pthread.h>

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
typedef std::stack<kgm_edge_descriptor> kgm_stack;
typedef std::vector<uint16_t> degree_stack;

ostream& operator<< (ostream& out, const kgm_adjacency_iterator& ai);
ostream& operator<< (ostream& out, const kgm_vertex_iterator& ai);
ostream& operator<< (ostream& out, const dfs_state& state);
ostream& operator<< (ostream& out, const std::pair<dfs_state,ugraph>& state);
ostream& operator<< (ostream& out, const dfs_stack& stack);

// Runtime

uint32_t KGM_GRAPH_SIZE; // shared
uint16_t KGM_UPPER_BOUND = 30; // shared
uint32_t KGM_LOWER_BOUND = 2; // shared, fixed
const uint32_t KGM_START_NODE = 0;
const uint64_t KGM_REPORT_INTERVAL = 0x10000000;
uint64_t KGM_REPORT_NEXT = KGM_REPORT_INTERVAL; // private (all thread reporting)
uint64_t KGM_STEPS = 0;

boost::scoped_ptr<boost::timer> KGM_TIMER; // shared, may be private for all
										   // threads to measure single thread span
bool running = true; // private

// TODO
int openMpThreadNumber();
// TODO
int openMpTotalThreads();

ugraph g; //shared
i_dfs_stack inactiveDfsStack; //private
dfs_stack dfsStack; //private
degree_stack degreeStack; //private


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

void sendFinish();
void broadcastNewSolution(int);

void dfs_step(
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
            std::cout << openMpNumber() << ": New spanning tree of degree: "<< degreeLimit << std::endl
                    << openMpNumber() << ": " << inactiveDfsStack << stack << std::endl
                    << openMpNumber() << ": of total size: " << inactiveDfsStack.size() + stack.size() << std::endl;

            if(degreeLimit == KGM_LOWER_BOUND)
                sendFinish();
            else
                broadcastNewSolution(degreeLimit);
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

void readInputFromFile(std::string filename) {
	std::ifstream in(filename.c_str());
	std::stringstream buffer;
	buffer << in.rdbuf();
	std::string contents(buffer.str());

	std::vector<std::string> lines;
	boost::split(lines, contents, boost::is_any_of("\n"));

	int m = lines.size()-2, n = lines[lines.size()-2].length();
	//std::cout << "loaded graph of size: " << m << "*" << n << std::endl;
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

bool hasExtraWork() {
    if (dfsStack.size() <= 1)
        return false;
    for (dfs_stack::iterator it = dfsStack.begin(); it < (dfsStack.end()-1); ++it)
    {
        if ((*it).v_it != (*it).v_it_end)
            return true;
    }
	return false;
}


// TODO Merge sendWork and GrabWork into single function grabWork, that operates
// on shared central (or thread-specific) stacks.
bool grabWork();

void sendWork(int processNumber) {
    int difference, diffRatio;
    kgm_vertex_iterator newVIt, newVIt_end;
    dfs_stack::iterator it;
    bool ok = false;

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
        return;
    (*it).v_it_end = newVIt;

    char* outputBuffer = new char [sizeof(uint16_t)*(KGM_GRAPH_SIZE*2 + 5)];
    uint16_t* outputBuffer16 = (uint16_t*) outputBuffer;
    int i = 1, cnt = 0;
    std::cout << MPI_MY_RANK << ": sending work - degree at level: " << it-dfsStack.begin()
            <<"( of"<< degreeStack.size() << ") = " << degreeStack.at(it-dfsStack.begin()) << std::endl;
    outputBuffer16[i++] = (uint16_t)(degreeStack.at(it-dfsStack.begin())); // spanning tree degree
    outputBuffer16[i++] = (uint16_t)(*newVIt); // starting state - v_it
    outputBuffer16[i++] = (uint16_t)(*newVIt_end); // starting state - v_it_end
    dfs_stack::iterator git;
    i_dfs_stack::iterator igit;
//    ++it;
    for (igit = inactiveDfsStack.begin(); igit != inactiveDfsStack.end(); ++igit)
    {
        outputBuffer16[i++] = (*igit).v;
        outputBuffer16[i++] = (*igit).a;
        ++cnt;
    }
    for (git = dfsStack.begin(); git != it; ++git)
    {
        outputBuffer16[i++] = (uint16_t)(*((*git).v_it));
        outputBuffer16[i++] = (uint16_t)(*((*git).a_it));
        ++cnt;
    }
    outputBuffer16[0] = (uint16_t)(cnt); // number of visited states - depth
    outputBuffer16[i++] = (uint16_t) 0; // terminating 00 bytes

    std::stringstream ss;
    ss << "[ ";
    for (int j = 0; j < i; ++j)
    {
        ss << outputBuffer16[j] << " ";
    }
    ss << "]";
    std::cout << MPI_MY_RANK << ": sent work - length (of uint16_t): " << i
            << " to " << processNumber << std::endl
            << ss.str() << std::endl;

    MPI_Send (outputBuffer, i*sizeof(uint16_t), MPI_CHAR, processNumber, MSG_WORK_SENT, MPI_COMM_WORLD);
    delete[] outputBuffer;

    /*std::cout << MPI_MY_RANK << ": my stack after send: " << std::endl
            << dfsStack << std::endl;*/
}

void acceptWork(MPI_Status& status) {
    int maxInputBufferSize = sizeof(uint16_t)*(KGM_GRAPH_SIZE*2 + 5), receivedNum;
    char* inputBuffer = new char [maxInputBufferSize];
    uint16_t* inputBuffer16 = (uint16_t*) inputBuffer;

    MPI_Recv(inputBuffer, maxInputBufferSize, MPI_CHAR, status.MPI_SOURCE, MSG_WORK_SENT, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, &receivedNum);

    int total16 = receivedNum / sizeof(uint16_t); // uint16_t index
    std::stringstream ss;
    ss << "[ ";
    for (int i = 0; i < total16; ++i)
    {
        ss << inputBuffer16[i] << " ";
    }
    ss << "]";
    std::cout << MPI_MY_RANK << ": received work - length (of uint16_t): " << total16
            << " from " << status.MPI_SOURCE << std::endl
            << MPI_MY_RANK << ": " << ss.str() << std::endl;

    if (inputBuffer16[total16-1] != 0)
    {
        std::cout << MPI_MY_RANK <<
                ": ERROR - acceptWork char* bad format - LAST UINT16_T != 0" << std::endl
                << "inputBuffer[total16-1] = " << inputBuffer[total16-1] << std::endl
                << "total16 = " << total16 << std::endl;
        throw("ERROR - acceptWork char* bad format - LAST UINT16_T != 0");
    }
    else {
        /*std::cout << MPI_MY_RANK <<
                ": OK - acceptWork char* valid format - LAST UINT16_T == 0" << std::endl;*/
    }

    int numberOfEdgesBefore = inputBuffer16[0];

    // Clear degree stack, put first state's degree into the stack
    degreeStack.clear();
    degreeStack.push_back(inputBuffer16[1]);

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

    // Initialize graph state from message
    int visitedEdgesCnt = 0;
//    g[inputBuffer16[2]].state = true;
    for (int j = 4, j_end = total16-1; j < j_end; j+=2)
    {
        i_dfs_state istate;
        istate.v = inputBuffer16[j];
        istate.a = inputBuffer16[j+1];
        inactiveDfsStack.push_back(istate);
        g[inputBuffer16[j]].state = true;
        g[inputBuffer16[j]].degree++;
        g[inputBuffer16[j+1]].state = true;
        g[inputBuffer16[j+1]].degree++;
        ++visitedEdgesCnt;
    }
    if (numberOfEdgesBefore != visitedEdgesCnt)
    {
        std::cout << MPI_MY_RANK <<
                ": ERROR - acceptWork char* bad format - inputBuffer16[0] != visited nodes" << std::endl
                << "numberOfEdgesBefore = inputBuffer16[0] = " << inputBuffer16[0] << std::endl
                << "visitedEdgesCnt = " << visitedEdgesCnt << std::endl;
        throw("ERROR - acceptWork char* bad format - inputBuffer16[0] != visited edges");
    }
    else
    {
        /*std::cout << MPI_MY_RANK <<
                ": OK - acceptWork char* valid format - inputBuffer16[0] == visited edges" << std::endl
                << "numberOfEdgesBefore = inputBuffer16[0] = " << inputBuffer16[0] << std::endl
                << "visitedEdgesCnt = " << visitedEdgesCnt << std::endl;*/
    }

    // Clear dfs stack and initialize first state
    dfsStack.clear();
    dfs_state firstState;
    tie(firstState.v_it, firstState.v_it_end) = vertices(g);
    firstState.v_it += inputBuffer16[2];
    firstState.v_it_end = firstState.v_it + (inputBuffer16[3] - inputBuffer16[2]);
    std::cout << MPI_MY_RANK << ": firstState.v_it_end = firstState.v_it + "
            << inputBuffer16[3]-inputBuffer16[2] << std::endl;
    tie(firstState.a_it, firstState.a_it_end) = adjacent_vertices(*(firstState.v_it),g);

    std::cout << MPI_MY_RANK << ": raw first state received: " << firstState << std::endl;

    if (!isValid_dfs_state(firstState, g))
    {
        std::cout << MPI_MY_RANK << ": is invalid" << std::endl;
        iterate_dfs_state(firstState, g);
//        std::cout << MPI_MY_RANK << ": firstState.a_it: " << (uint16_t)(*(firstState.a_it)) << std::endl;
//        std::cout << MPI_MY_RANK << ": firstState.a_it_end - a_it: " << (uint16_t)(firstState.a_it_end - firstState.a_it) << std::endl;
        std::cout << MPI_MY_RANK << ": iterated to: " << firstState << std::endl;
    } else {
        std::cout << MPI_MY_RANK << ": is valid" << std::endl;
    }
    if (!isValid_dfs_state(firstState, g))
    {
        std::cout << MPI_MY_RANK << ": is invalid after iteration" << std::endl
                  << MPI_MY_RANK << ": has accepted no valid work!" << std::endl;
    }
    else
    {
        std::cout << MPI_MY_RANK << ": successfully accepted work" << std::endl
                << "   dfs stack size: " << dfsStack.size() << std::endl
                << "   degree: " << degreeStack.back() << std::endl
                << "   edges before this state: " << visitedEdgesCnt << std::endl;
        dfsStack.push_back(firstState);
    }

}


void divideWork() {
    std::cout << ": divideWork()" << std::endl;
	int processNumber = 1; // starting with the next process
	while(hasExtraWork() && processNumber < openMpTotalThreads()) {
		sendWork(processNumber);
		++processNumber;
	}
}

void initStack() {
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


void iterateStack() {

	while(running) {
		while (!dfsStack.empty())
		{
			++KGM_STEPS;

			if (dfsStack.empty())
				break;

			dfs_step(dfsStack, degreeStack, g, KGM_UPPER_BOUND);

			if (KGM_STEPS >= KGM_REPORT_NEXT)
			{
				std::cout << ": Stack report at " << KGM_TIMER->elapsed() << ":" << std::endl
						<< inactiveDfsStack << dfsStack << std::endl;
				KGM_REPORT_NEXT += KGM_REPORT_INTERVAL;
			}
		}

		if (!grabWork())
			running = false;
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

	readInputFromFile(filename);

    initStack();


    KGM_TIMER.reset(new boost::timer);

    iterateStack();

	std::cout << ": ***** ENDED in time " << KGM_TIMER->elapsed() << std::endl;

	return 0;
}

