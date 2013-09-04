#ifndef __SOINN_H__
#define __SOINN_H__

#include <fstream>
#include <iostream>
#include <string>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adj_list_serialize.hpp>
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
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version)
		{
			ar & BOOST_SERIALIZATION_NVP(weight);
			ar & BOOST_SERIALIZATION_NVP(class_id);
			ar & BOOST_SERIALIZATION_NVP(M);
			ar & BOOST_SERIALIZATION_NVP(error);
			ar & BOOST_SERIALIZATION_NVP(R);
		}
	};

	struct EdgeProperties
	{
		int age;
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version)
		{
			ar & BOOST_SERIALIZATION_NVP(age);
		}
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
		void save(std::string filename);
		void load(std::string filename);
		void clear();
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
	private:
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version)
		{
			ar & BOOST_SERIALIZATION_NVP(class_count);
			ar & BOOST_SERIALIZATION_NVP(iteration_count);
			ar & BOOST_SERIALIZATION_NVP(age_max);
			ar & BOOST_SERIALIZATION_NVP(lambda);
			ar & BOOST_SERIALIZATION_NVP(alpha1); 
			ar & BOOST_SERIALIZATION_NVP(alpha2); 
			ar & BOOST_SERIALIZATION_NVP(alpha3);
			ar & BOOST_SERIALIZATION_NVP(beta); 
			ar & BOOST_SERIALIZATION_NVP(gamma);
			ar & BOOST_SERIALIZATION_NVP(C);
			ar & BOOST_SERIALIZATION_NVP(graph);
		}
	};
}

#endif /*SOINN.h*/