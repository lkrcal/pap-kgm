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
#include "mpi.h"

using namespace std;
//using namespace boost::filesystem;
using namespace boost;

static int32_t MPI_MY_RANK; // my process number
static int32_t MPI_PROCESSES; // number of processes
static int32_t MPI_REQUEST_PROCESS; // process from which we will be requesting work

enum state { // states of the process
	SLEEPING,
	WORKING,
	NEED_WORK,
	FINISHED
};

#define CHECK_MSG_AMOUNT  100

#define MSG_WORK_REQUEST 1000
#define MSG_WORK_SENT    1001
#define MSG_WORK_NOWORK  1002
#define MSG_TOKEN        1003
#define MSG_FINISH       1004
#define MSG_NEW_SOLUTION	1005


static int32_t KGM_GRAPH_SIZE;
static int32_t KGM_UPPER_BOUND = 30;
static int32_t KGM_LOWER_BOUND = 2;
static int32_t KGM_START_NODE = 0;
static uint64_t KGM_REPORT_INTERVAL = 0x10000000;
static state PROCESS_STATE = SLEEPING;
static bool running = true;
//static path graphSource;

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
typedef std::vector<dfs_state> dfs_stack;
typedef std::stack<kgm_edge_descriptor> kgm_stack;
typedef std::vector<uint16_t> degree_stack;

ugraph g;

ostream& operator<< (ostream& out, const kgm_adjacency_iterator& ai);
ostream& operator<< (ostream& out, const kgm_vertex_iterator& ai);
ostream& operator<< (ostream& out, const dfs_state& state);
ostream& operator<< (ostream& out, const std::pair<dfs_state,ugraph>& state);
ostream& operator<< (ostream& out, const dfs_stack& stack);

bool iterate_dfs_state(dfs_state& dfsState, const ugraph& graph)
{
    for (++dfsState.a_it; dfsState.a_it != dfsState.a_it_end; ++dfsState.a_it)
    {
        if (graph[*(dfsState.a_it)].state == false)
            return true;
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

void dfs_step(
        dfs_stack& stack,
        degree_stack& dstack,
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

        dstack.push_back(std::max(graph[*(stack.back().v_it)].degree, dstack.back()));

        if (dstack.back() >= degreeLimit)
        {   // not going to be a better solution
            dstack.pop_back();
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
            std::cout << "New spanning tree - degree: "<< dstack.back() << std::endl;
            degreeLimit = dstack.back();
            // TODO possible spanning tree improvement
            std::cout << stack << std::endl;
        }
        return;
    }

    dfs_state prev (stack.back());
    graph[*(prev.v_it)].degree -= 1;
    graph[*(prev.a_it)].degree = 0;
    graph[*(prev.a_it)].state = false;

    dstack.pop_back();

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

bool hasExtraWork() {
	// TODO - asi Lubos? :o) detekce toho, jestli mame co delit
	return false;
}

void sendWork() {
	// TODO - asi Lubos? :o) serializovat a poslat
}

void acceptWork(char* _buffer, int _buffer_size) {
	// TODO - asi Lubos? :o) deserializovat a prijmout
}

void requestWork() {
	int message = 1;
	std::cout << MPI_MY_RANK << " is requesting work from " << MPI_REQUEST_PROCESS;

	MPI_Send(&message, 1, MPI_INT, MPI_REQUEST_PROCESS, MSG_WORK_REQUEST, MPI_COMM_WORLD);

	++MPI_REQUEST_PROCESS % MPI_PROCESSES;
	if(MPI_REQUEST_PROCESS == MPI_MY_RANK) {
		MPI_REQUEST_PROCESS++;
	}
}

void receiveMessage() {
	int flag = 0;
	MPI_Status status;

	MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
	if (flag)
	{
	  int buffer_size;
	  char* buffer = new char[buffer_size];
	  MPI_Recv(buffer, buffer_size, MPI_CHAR, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	  switch (status.MPI_TAG)
	  {
		 case MSG_WORK_REQUEST :
			 MPI_Recv(&buffer_size, 1, MPI_INT,  MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			 if(hasExtraWork()) {
				 sendWork();
			 } else {
				 // not enough work to send, sending refusal
				 MPI_Send (buffer, 0, MPI_INT, status.MPI_SOURCE, MSG_WORK_NOWORK, MPI_COMM_WORLD);
			 }
			 break;
		 case MSG_WORK_SENT :
			  acceptWork(buffer, buffer_size);
			  break;
		 case MSG_WORK_NOWORK :
			  break;
		 case MSG_TOKEN :
			  break;
		 case MSG_NEW_SOLUTION :
			  break;
		 case MSG_FINISH :
		      MPI_Finalize();
		      exit (0);
		      break;
		 default : // error
			 break;
	  }
	}
}

void work() {
	while(running) {
		switch(PROCESS_STATE) {
			case WORKING:
				break;
			case NEED_WORK:
				break;
			case FINISHED:
				break;
			default:
				throw "Unknown state detected";
		}
	}
}

void printResult(double time, uint64_t steps) {
	std::cout << "-------------" << std::endl;
	std::cout << "TOTAL STEPS: " << steps << std::endl;
	std::cout << "In time: " << time << std::endl;
}

void iterateStack() {
	dfs_state firstState;
	g[KGM_START_NODE].state = true;
	if (!create_dfs_state(firstState,g))
	{
		std::cerr << "Failed to initialize first state" << std::endl;
		exit(-4);
	}

	dfs_stack stack;
	stack.push_back(firstState);
	degree_stack dstack; // max degree
	dstack.push_back(0);

	uint16_t minSPdegree = KGM_UPPER_BOUND;
	uint64_t steps = 0;
	uint64_t psteps = KGM_REPORT_INTERVAL;
	boost::scoped_ptr<boost::timer> timer (new boost::timer);

	while (!stack.empty())
	{
		if ((psteps % CHECK_MSG_AMOUNT)==0)
		{
			receiveMessage();
		}

		dfs_step(stack, dstack, g, minSPdegree);
		if (steps >= psteps)
		{
			std::cout << "Stack at " << timer->elapsed() << ":" << std::endl
					<< make_pair(stack.front(),g) << std::endl;
			psteps += KGM_REPORT_INTERVAL;
		}
		++steps;
	}
	double time = timer->elapsed();
	timer.reset();

	printResult(time, steps);
}

int main(int argc, char ** argv) {

    if (argc <= 1)
    {
        std::cerr << "Not enough arguments." << std::endl;
        return -1;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_MY_RANK); // my rank
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_PROCESSES); // number of processes
    MPI_Barrier(MPI_COMM_WORLD); // loading all processes

    std::string filename (argv[1]);
	if (filename.empty())
		return -1;

	readInputFromFile(filename);

    if(MPI_MY_RANK == 0) {
		iterateStack();
    }

    MPI_Finalize();

    return 0;
}
