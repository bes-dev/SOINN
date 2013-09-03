#ifndef __SOINN_H__
#define __SOINN_H__

#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace soinn
{
	struct VertexProperties
	{
		boost::numeric::ublas::vector<double> weight;
		int class_id;
		double M;
		double error;
		double R;
	};

	struct EdgeProperties
	{
		int age;
	};

	typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, boost::property<boost::vertex_index_t, size_t, VertexProperties>, EdgeProperties> Graph;
	typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef boost::graph_traits<Graph>::edge_descriptor Edge;
	typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
	typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;
	typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIterator;
	typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIterator;
	typedef boost::graph_traits<Graph>::vertices_size_type VerticesSizeType;
	typedef std::map<Vertex, int> ComponentMap;
	
	class SOINN
	{
	public:
		SOINN(int lambda = 30, int age_max = 10, double C = 1., double alpha1 = 0.16, double alpha2 = 0.25, double alpha3 = 0.25, double beta = 0.67, double gamma = 0.75);
		~SOINN();

		void init(boost::numeric::ublas::vector<double> v1, boost::numeric::ublas::vector<double> v2);
		int addSignal(const boost::numeric::ublas::vector<double> &x);
		void classify();
		Graph getGraph();
	private:
		Graph graph;
		int class_count;
		int iteration_count;
		int age_max;
		int lambda;
		double alpha1, alpha2, alpha3;
		double beta, gamma;
		double C;

		int findBestMatches(const boost::numeric::ublas::vector<double> &x, VertexIterator &match_first, VertexIterator &match_second);
		int incrementEdgeAge(const VertexIterator &vertex);
		void updateWeights(const boost::numeric::ublas::vector<double> &x, const VertexIterator &vertex);
		void deleteOldEdges();
		void addNewNodeAndRemoveUnnecessaryNodes();
		Vertex findMaxErrorNode();
		Vertex findMaxLocErrorNode(const Vertex &vertex);
		double distance(const boost::numeric::ublas::vector<double> &a, const boost::numeric::ublas::vector<double> &b);
		double getSimilarityThreshold(const VertexIterator &vertex);
		double getMeanM();
		bool isWithinThreshold(const boost::numeric::ublas::vector<double> &x, const VertexIterator &match_first, const VertexIterator &match_second);
	};
}

#endif /*SOINN.h*/