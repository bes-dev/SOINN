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

SOINN::SOINN(int lambda, int age_max, double C, double alpha1, double alpha2, double alpha3, double beta, double gamma): class_count(0), iteration_count(0) ,lambda(lambda), age_max(age_max), alpha1(alpha1), alpha2(alpha2), alpha3(alpha3), beta(beta), gamma(gamma)
{
}

SOINN::~SOINN()
{
}

void SOINN::init(boost::numeric::ublas::vector<double> v1, boost::numeric::ublas::vector<double> v2)
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

int SOINN::addSignal(const boost::numeric::ublas::vector<double> &x)
{
	iteration_count++;
	VertexIterator match_first, match_second;
	findBestMatches(x, match_first, match_second);
	if(!isWithinThreshold(x, match_first, match_second))
	{
		Vertex vertex = boost::add_vertex(graph);
		graph[vertex].weight = boost::numeric::ublas::vector<double>(x);
		graph[vertex].M = 1.;
		graph[vertex].R = 0.;
		graph[vertex].error = 0.;
		graph[vertex].class_id = -1;
		return 0;
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
	incrementEdgeAge(match_first);
	graph[*match_first].error += distance(x, graph[*match_first].weight);
	graph[*match_first].M++;
	updateWeights(x, match_first);
	deleteOldEdges();
	if(iteration_count % lambda == 0)
	{
		addNewNodeAndRemoveUnnecessaryNodes();
	}

	return 0;
}

double SOINN::distance(const boost::numeric::ublas::vector<double> &a, const boost::numeric::ublas::vector<double> &b)
{
	return ublas::norm_2(a - b);
}

int SOINN::findBestMatches(const boost::numeric::ublas::vector<double> &x, VertexIterator &match_first, VertexIterator &match_second)
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

double SOINN::getSimilarityThreshold(const VertexIterator &vertex)
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

bool SOINN::isWithinThreshold(const boost::numeric::ublas::vector<double> &x, const VertexIterator &match_first, const VertexIterator &match_second)
{
	if(distance(x, graph[*match_first].weight) > getSimilarityThreshold(match_first))
	{
		return false;
	}
	if(distance(x, graph[*match_second].weight) > getSimilarityThreshold(match_second))
	{
		return false;
	}
	return true;
}

int SOINN::incrementEdgeAge(const VertexIterator &vertex)
{
	OutEdgeIterator edge_out_current, edge_out_end;
	boost::tie(edge_out_current, edge_out_end) = boost::out_edges(*vertex, graph);
	for(; edge_out_current != edge_out_end; edge_out_current++)
	{
		graph[*edge_out_current].age++;
	}
	return 0;
}

void SOINN::updateWeights(const boost::numeric::ublas::vector<double> &x, const VertexIterator &vertex)
{
	graph[*vertex].weight += E1(graph[*vertex].M) * (x - graph[*vertex].weight);
	AdjacencyIterator vertex_current, vertex_end;
	boost::tie(vertex_current, vertex_end) = boost::adjacent_vertices(*vertex, graph);
	for(; vertex_current != vertex_end; vertex_current++)
	{
		graph[*vertex_current].weight += E2(graph[*vertex_current].M) * (x - graph[*vertex_current].weight);
	}
}

void SOINN::deleteOldEdges()
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

void SOINN::addNewNodeAndRemoveUnnecessaryNodes()
{
	Vertex q = findMaxErrorNode();
	Vertex f = findMaxLocErrorNode(q);
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
	double mean_m = getMeanM();
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

Vertex SOINN::findMaxErrorNode()
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

Vertex SOINN::findMaxLocErrorNode(const Vertex &vertex)
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

double SOINN::getMeanM()
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

void SOINN::classify()
{
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

Graph SOINN::getGraph()
{
	return graph;
}