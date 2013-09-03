#include <iostream>
#include <vector>
#include <string>
#include <boost/numeric/ublas/vector.hpp>
#include "SOINN.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <random>
#include <ctime>
#include <boost/foreach.hpp>
#include <boost/utility.hpp> 
#include <boost/graph/connected_components.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/adjacency_list.hpp>

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
	//std::srand(0);
	soinn::SOINN model;
	ublas::vector<double> a(2);
	ublas::vector<double> x1(2), x2(2);
	x1(0) = 320;
	x1(1) = 220;
	x2(0) = 320;
	x2(1) = 240;

	model.init(x1, x2);
	int sign = 1;
	for(int i = 0; i < 4000; i++)
	{
		std::cout<<i<<"\n";

		a(0) = 200 + std::rand()%200;
		a(1) = 200 + std::rand()%200;

		//int rnd = std::rand()%3;
		//if(rnd == 0)
		//{
		//	a(0) = 300 + std::rand()%40;
		//	a(1) = 220 + std::rand()%40;
		//}
		//else if(rnd == 1)
		//{
		//	a(0) = 420 + std::rand()%40;
		//	a(1) = 290 + std::rand()%40;
		//}
		//else
		//{
		//	a(0) = 100 + std::rand()%40;
		//	a(1) = 120 + sign*std::sqrt(400-(a(0)-120)*(a(0)-120));
		//	sign *= -1;
		//}

		model.addSignal(a);
	}
	model.classify();
	cv::Mat img(480, 640, CV_32FC3);
	cv::rectangle(img, cv::Point(200, 200), cv::Point(400, 400), cv::Scalar(0, 255, 0));
	//cv::rectangle(img, cv::Point(300, 220), cv::Point(340, 260), cv::Scalar(0, 255, 0));
	//cv::rectangle(img, cv::Point(420, 290), cv::Point(460, 330), cv::Scalar(0, 255, 0));
	//cv::circle(img, cv::Point(120,120), 20, cv::Scalar(0, 255, 0));
	draw(img, model.getGraph());
	cv::imshow("img", img);
	cv::waitKey();
	return 0;
}