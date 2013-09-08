/*
* SOINN.cpp
*
*      Author: Sergei Belousov aka BeS
*/

#include "SOINN.h"

#include <boost/foreach.hpp>
#include <boost/utility.hpp> 
#include <boost/graph/connected_components.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/iteration_macros.hpp>

using namespace soinn;
using namespace boost::numeric;

#define E1(t) 1./t
#define E2(t) 1./(100*t)

SOINN::SOINN(int lambda, int age_max, double C, double alpha1, double alpha2, double alpha3, double beta, double gamma): layer_flag(1), class_count(0), iteration_count(0) ,lambda(lambda), age_max(age_max), C(C), alpha1(alpha1), alpha2(alpha2), alpha3(alpha3), beta(beta), gamma(gamma)
{
	std::srand(std::time(0));
}

SOINN::~SOINN()
{
}

void SOINN::init(boost::numeric::ublas::vector<double> v1, boost::numeric::ublas::vector<double> v2)
{
	initGraph(first_layer, v1, v2);
}

void SOINN::initGraph(Graph &graph, boost::numeric::ublas::vector<double> v1, boost::numeric::ublas::vector<double> v2)
{
	Vertex vertex_v1 = boost::add_vertex(graph);
	graph[vertex_v1].weight = boost::numeric::ublas::vector<double>(v1);
	graph[vertex_v1].M = 1.;
	graph[vertex_v1].R = 0.;
	graph[vertex_v1].error = 0.;
	graph[vertex_v1].class_id = -1;
	Vertex vertex_v2 = boost::add_vertex(graph);
	graph[vertex_v2].weight = boost::numeric::ublas::vector<double>(v2);
	graph[vertex_v2].M = 1.;
	graph[vertex_v2].R = 0.;
	graph[vertex_v2].error = 0.;
	graph[vertex_v2].class_id = -1;
}

void SOINN::addSignal(const boost::numeric::ublas::vector<double> &x)
{
	addSignalInGraph(first_layer, x);
}

void SOINN::addSignalInGraph(Graph &graph, const boost::numeric::ublas::vector<double> &x)
{
	iteration_count++;
	VertexIterator match_first, match_second;
	findBestMatches(graph, x, match_first, match_second);
	if(!isWithinThreshold(graph, x, match_first, match_second))
	{
		Vertex vertex = boost::add_vertex(graph);
		graph[vertex].weight = boost::numeric::ublas::vector<double>(x);
		graph[vertex].M = 1.;
		graph[vertex].R = 0.;
		graph[vertex].error = 0.;
		graph[vertex].class_id = -1;
		return;
	}
	std::pair<Edge, bool> pe = boost::edge(*match_first, *match_second, graph);
	if(!pe.second)
	{
		Edge e = boost::add_edge(*match_first, *match_second, graph).first;
		graph[e].age = 0;
	}
	else
	{
		graph[pe.first].age = 0;
	}
	incrementEdgeAge(graph, match_first);
	graph[*match_first].error += distance(x, graph[*match_first].weight);
	graph[*match_first].M++;
	updateWeights(graph, x, match_first);
	deleteOldEdges(graph);
	if(iteration_count % lambda == 0)
	{
		addNewNodeAndRemoveUnnecessaryNodes(graph);
	}
}

double SOINN::distance(const boost::numeric::ublas::vector<double> &a, const boost::numeric::ublas::vector<double> &b)
{
	return ublas::norm_2(a - b);
}

int SOINN::findBestMatches(Graph &graph, const boost::numeric::ublas::vector<double> &x, VertexIterator &match_first, VertexIterator &match_second)
{
	double dist_first = std::numeric_limits<double>::max();
	double dist_second = std::numeric_limits<double>::max();

	VertexIterator vertex_current, vertex_end;
	boost::tie(vertex_current, vertex_end) = boost::vertices(graph);

	for(; vertex_current != vertex_end; ++vertex_current)
	{
		double dist = distance(x, graph[*vertex_current].weight);
		if(dist < dist_first)
		{
			match_second = match_first;
			dist_second = dist_first;
			match_first = vertex_current;
			dist_first = dist;
		}
		else if(dist < dist_second)
		{
			match_second = vertex_current;
			dist_second = dist;
		}
	}

	return 0;
}

double SOINN::getSimilarityThreshold(Graph &graph, const VertexIterator &vertex)
{
	double dist = 0.0;
	if(boost::out_degree(*vertex, graph) == 0)
	{
		dist = std::numeric_limits<double>::max();
		VertexIterator vertex_current, vertex_end;
		boost::tie(vertex_current, vertex_end) = boost::vertices(graph);
		for(; vertex_current != vertex_end; ++vertex_current)
		{
			if(vertex_current != vertex)
			{
				double dist_current = distance(graph[*vertex].weight, graph[*vertex_current].weight);
				if(dist_current < dist)
				{
					dist = dist_current;
				}
			}
		}		
	}
	else
	{
		dist = std::numeric_limits<double>::min();
		AdjacencyIterator vertex_current, vertex_end;
		boost::tie(vertex_current, vertex_end) = boost::adjacent_vertices(*vertex, graph);
		for(; vertex_current != vertex_end; ++vertex_current)
		{
			double dist_current = distance(graph[*vertex].weight, graph[*vertex_current].weight);
			if(dist_current > dist)
			{
				dist = dist_current;
			}
		}
	}

	return dist;
}

bool SOINN::isWithinThreshold(Graph &graph, const boost::numeric::ublas::vector<double> &x, const VertexIterator &match_first, const VertexIterator &match_second)
{
	double T1 = layer_flag ? getSimilarityThreshold(graph, match_first) : Tc;
	double T2 = layer_flag ? getSimilarityThreshold(graph, match_second) : Tc;
	if(distance(x, graph[*match_first].weight) > T1)
	{
		return false;
	}
	if(distance(x, graph[*match_second].weight) > T2)
	{
		return false;
	}
	return true;
}

int SOINN::incrementEdgeAge(Graph &graph, const VertexIterator &vertex)
{
	OutEdgeIterator edge_out_current, edge_out_end;
	boost::tie(edge_out_current, edge_out_end) = boost::out_edges(*vertex, graph);
	for(; edge_out_current != edge_out_end; edge_out_current++)
	{
		graph[*edge_out_current].age++;
	}
	return 0;
}

void SOINN::updateWeights(Graph &graph, const boost::numeric::ublas::vector<double> &x, const VertexIterator &vertex)
{
	graph[*vertex].weight += E1(graph[*vertex].M) * (x - graph[*vertex].weight);
	AdjacencyIterator vertex_current, vertex_end;
	boost::tie(vertex_current, vertex_end) = boost::adjacent_vertices(*vertex, graph);
	for(; vertex_current != vertex_end; vertex_current++)
	{
		graph[*vertex_current].weight += E2(graph[*vertex_current].M) * (x - graph[*vertex_current].weight);
	}
}

void SOINN::deleteOldEdges(Graph &graph)
{
	EdgeIterator edge_current, edge_end;
	boost::tie(edge_current, edge_end) = boost::edges(graph);
	EdgeIterator edge_next = edge_current;
	for(; edge_current != edge_end; edge_current = edge_next)
	{
		edge_next ++;
		if(graph[*edge_current].age > age_max)
		{
			Vertex vertex_s = boost::source(*edge_current, graph);
			Vertex vertex_t = boost::target(*edge_current, graph);
			boost::remove_edge(vertex_s, vertex_t, graph);
		}
	}	
}

void SOINN::addNewNodeAndRemoveUnnecessaryNodes(Graph &graph)
{
	Vertex q = findMaxErrorNode(graph);
	Vertex f = findMaxLocErrorNode(graph ,q);
	Vertex vertex = boost::add_vertex(graph);
	graph[vertex].weight = (graph[q].weight + graph[f].weight) * 0.5;
	graph[vertex].error = alpha1 * (graph[q].error + graph[f].error);
	graph[vertex].M = alpha2 * (graph[q].M + graph[f].M);
	double Rq = graph[q].error / graph[q].M;
	double Rf = graph[f].error / graph[f].M;
	graph[vertex].R = alpha3 * ( Rq + Rf );
	double Eq = beta * graph[q].error;
	double Ef = beta * graph[f].error;
	double Mq = gamma * graph[q].M;
	double Mf = gamma * graph[f].M;
	if(Eq / Mq < Rq && Ef / Mf < Rf)
	{
		graph[q].error = Eq;
		graph[q].M = Mq;
		graph[q].R = Rq;
		graph[f].error = Ef;
		graph[f].M = Mf;
		graph[f].R = Rf;
		Edge vq = boost::add_edge(q, vertex, graph).first;
		Edge vf = boost::add_edge(vertex, f, graph).first;
		graph[vq].age = 0;
		graph[vf].age = 0;
		boost::remove_edge(q, f, graph);
	}
	else
	{
		boost::clear_vertex(vertex, graph);
		boost::remove_vertex(vertex, graph);
	}
	double mean_m = getMeanM(graph);
	VertexIterator vertex_current, vertex_end;
	boost::tie(vertex_current, vertex_end) = boost::vertices(graph);
	VertexIterator vertex_next;
	for(vertex_next = vertex_current; vertex_current != vertex_end; vertex_current = vertex_next)
	{
		vertex_next++;
		if(boost::num_vertices(graph) == 2)
		{
			return;
		}
		if((boost::out_degree(*vertex_current, graph) == 1 && graph[*vertex_current].M < C * mean_m) || (boost::out_degree(*vertex_current, graph) == 0))
		{
			boost::clear_vertex(*vertex_current, graph);
			boost::remove_vertex(*vertex_current, graph);
		}
	}
}

Vertex SOINN::findMaxErrorNode(Graph &graph)
{
	Vertex vertex_max_error;
	double max_error = 0;
	VertexIterator vertex_current, vertex_end;
	boost::tie(vertex_current, vertex_end) = boost::vertices(graph);
	for(; vertex_current != vertex_end; vertex_current++)
	{
		if(graph[*vertex_current].error > max_error)
		{
			vertex_max_error = *vertex_current;
			max_error = graph[*vertex_current].error;
		}
	}
	return vertex_max_error;
}

Vertex SOINN::findMaxLocErrorNode(Graph &graph, const Vertex &vertex)
{
	Vertex vertex_max_error = vertex;
	double max_error = 0;
	AdjacencyIterator vertex_current, vertex_end;
	boost::tie(vertex_current, vertex_end) = boost::adjacent_vertices(vertex, graph);
	for(; vertex_current != vertex_end; ++vertex_current)
	{
		if(graph[*vertex_current].error > max_error)
		{
			vertex_max_error = *vertex_current;
			max_error = graph[*vertex_current].error;
		}
	}
	return vertex_max_error;
}

double SOINN::getMeanM(Graph &graph)
{
	double mean = 0;
	VertexIterator vertex_current, vertex_end;
	boost::tie(vertex_current, vertex_end) = boost::vertices(graph);
	for(; vertex_current != vertex_end; ++vertex_current)
	{
		mean += graph[*vertex_current].M;
	}
	mean /= boost::num_vertices(graph);
	return mean;
}

void SOINN::classifyGraph(Graph &graph)
{
	addNewNodeAndRemoveUnnecessaryNodes(graph);
	size_t index = 0;
	BGL_FORALL_VERTICES(v, graph, Graph)
	{
		boost::put(boost::vertex_index, graph, v, index++);
	}
	ComponentMap component;
	boost::associative_property_map<ComponentMap> component_map(component);
	class_count = connected_components(graph, component_map);
	BGL_FORALL_VERTICES(v, graph, Graph)
	{
		graph[v].class_id = boost::get(component_map, v);
	}
}

void SOINN::classify(int iterations)
{
	addNewNodeAndRemoveUnnecessaryNodes(first_layer);
	classifyGraph(first_layer);
	iteration_count = 0;
	trainSecondLayer(iterations);
	iteration_count = 0;
	classifyGraph(second_layer);
}

Graph SOINN::getFirstLayer()
{
	return first_layer;
}

Graph SOINN::getSecondLayer()
{
	return second_layer;
}

void SOINN::clearGraph(Graph &graph)
{
	graph.clear();
}

void SOINN::clear()
{
	clearGraph(first_layer);
}

void SOINN::save(std::string filename)
{
   	std::ofstream ofs(filename.c_str());
	boost::archive::xml_oarchive oa(ofs);
	oa << BOOST_SERIALIZATION_NVP(this);
}

void SOINN::load(std::string filename)
{
	clear();
   	std::ifstream ifs(filename.c_str());
	boost::archive::xml_iarchive ia(ifs);
	ia >> BOOST_SERIALIZATION_NVP(*const_cast<SOINN*>(this));
}

void SOINN::calcT()
{
	double dw = 0.;
	EdgeIterator edge_begin, edge_end;
	boost::tie(edge_begin, edge_end) = boost::edges(first_layer);
	for(EdgeIterator i = edge_begin; i != edge_end; ++i)
	{
		dw += distance(first_layer[i->m_source].weight, first_layer[i->m_target].weight);
	}
	dw /= boost::num_edges(first_layer);
	std::map<int, std::map<int, double>> db;
	for(int i = 0; i < class_count; ++i)
	{
		for(int j = 0; j < class_count; ++j)
		{
			db[i][j] = 0.;
		}
	}
	VertexIterator vertex_begin, vertex_end;
	boost::tie(vertex_begin, vertex_end) = boost::vertices(first_layer);
	for(VertexIterator i = vertex_begin; i != vertex_end; ++i)
	{
		for(VertexIterator j = vertex_begin; j != vertex_end; ++j)
		{
			if(first_layer[*i].class_id != first_layer[*j].class_id)
			{
				double dist_tmp = distance(first_layer[*i].weight, first_layer[*j].weight);
				if(dist_tmp > db[first_layer[*i].class_id][first_layer[*j].class_id])
				{
					db[first_layer[*i].class_id][first_layer[*j].class_id] = dist_tmp;
				}
			}
		}
	}
	Tc = std::numeric_limits<double>::max();
	for(int i = 0; i < class_count; ++i)
	{
		for(int j = 0; j < class_count; ++j)
		{
			if(db[i][j] > dw && db[i][j] < Tc)
			{
				Tc = db[i][j];
			}
		}
	}
}

void SOINN::trainSecondLayer(int iterations)
{
	layer_flag = 0;

	clearGraph(second_layer);
	calcT();
	std::vector<VertexIterator> data;
	VertexIterator vertex_begin, vertex_end;
	boost::tie(vertex_begin, vertex_end) = boost::vertices(first_layer);
	for(VertexIterator i = vertex_begin; i != vertex_end; ++i)
	{
		data.push_back(i);
	}
	initGraph(second_layer, first_layer[*data[0]].weight, first_layer[*data[1]].weight);
	std::random_shuffle(data.begin(), data.end());
	int size = data.size();
	for(int i = 0; i < iterations; ++i)
	{
		std::cout<<i<<"\n";
		int id = std::rand() % size;
		addSignalInGraph(second_layer, first_layer[*data[id]].weight);
	}

	layer_flag = 1;
}