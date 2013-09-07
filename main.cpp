#include <iostream>
#include <vector>
#include <string>
#include <boost/numeric/ublas/vector.hpp>
#include "SOINN.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <random>
#include <ctime>
#include <fstream>
#include <iostream>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/foreach.hpp>
#include <boost/utility.hpp> 
#include <boost/graph/connected_components.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adj_list_serialize.hpp>

using namespace boost::numeric;
using namespace soinn;

void draw(cv::Mat &image, Graph graph)
{
	const cv::Scalar colors[6] = {cv::Scalar(255, 0, 0), cv::Scalar(0, 255, 0), cv::Scalar(0, 0, 255), cv::Scalar(255, 255, 0), cv::Scalar(0, 255, 255), cv::Scalar(255, 0, 255)};
	EdgeIterator edge_begin, edge_end;
	boost::tie(edge_begin, edge_end) = boost::edges(graph);
	for(EdgeIterator i = edge_begin; i != edge_end; i++)
	{
		Vertex vertex_s = boost::source(*i, graph);
		Vertex vertex_t = boost::target(*i, graph);
		cv::Point p1(graph[vertex_s].weight[0], graph[vertex_s].weight[1]);
		cv::Point p2(graph[vertex_t].weight[0], graph[vertex_t].weight[1]);
		cv::line(image, p1, p2, colors[graph[vertex_s].class_id%6]);
	}
	VertexIterator vertex_begin, vertex_end;
	boost::tie(vertex_begin, vertex_end) = boost::vertices(graph);
	for(VertexIterator i = vertex_begin; i != vertex_end; i++)
	{
		cv::circle(image, cv::Point(graph[*i].weight[0], graph[*i].weight[1]), 3, colors[graph[*i].class_id%6]);
	}
}


int main()
{
	std::srand(std::time(0));
	soinn::SOINN model;
	ublas::vector<double> a(2);
	ublas::vector<double> x1(2), x2(2);
	x1(0) = 320;
	x1(1) = 220;
	x2(0) = 320;
	x2(1) = 240;

	model.init(x1, x2);
	cv::Mat img = cv::imread("img.png", 0);
	std::vector<cv::Point> points;
	for(int i = 0; i < img.size().width; ++i)
	{
		for(int j = 0; j < img.size().height; ++j)
		{
			if(img.at<char>(cv::Point(i, j)) != 0)
			{
				points.push_back(cv::Point(i, j));
			}
		}
	}

	int size = points.size();
	for(int i = 0; i < 10000; ++i)
	{
		std::cout<<i<<"\n";
		int id = std::rand() % size;
		a(0) = points[id].x;
		a(1) = points[id].y;
		model.addSignal(a);
	}

	model.classify();
	//soinn::SOINN model;
	//model.load("filename.xml");

	cv::imshow("img", img);

	cv::Mat img1(480, 640, CV_32FC3);
	//cv::rectangle(img1, cv::Point(200, 200), cv::Point(400, 400), cv::Scalar(0, 255, 0));
	//cv::rectangle(img1, cv::Point(300, 220), cv::Point(340, 260), cv::Scalar(0, 255, 0));
	//cv::rectangle(img1, cv::Point(420, 290), cv::Point(460, 330), cv::Scalar(0, 255, 0));
	//cv::circle(img1, cv::Point(120,120), 20, cv::Scalar(0, 255, 0));

	draw(img1, model.getFirstLayer());
	//draw(img1, model.getGraph());
	cv::imshow("img1", img1);


	cv::Mat img2(480, 640, CV_32FC3);
	draw(img2, model.getSecondLayer());
	cv::imshow("img2", img2);

//	std::cout<<boost::num_vertices(model.getFirstLayer())<<"; "<<boost::num_vertices(model.getSecondLayer());

	model.save("filename.xml");

	cv::waitKey();
	return 0;
}