#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#include<mpi.h>
#include <map>
#include <map>
#include <map>
#pragma once

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>

using namespace std;
using namespace msclr::interop;

int IMAGE_BYTES, CLUSTER_BYTES;
int* bytes_color, * ByteCentroid;
int nCentroids, nIterations, Byte_size;
const char* inputFile, * outputFile;

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* input;


	int OriginalImageWidth, OriginalImageHeight;

	//*********************************************************Read Image and save it to local arrayss*************************	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int* Red = new int[BM.Height * BM.Width];
	int* Green = new int[BM.Height * BM.Width];
	int* Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height * BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i * BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values
		}

	}
	return input;
}
void createImage(int* image, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);
	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			if (image[i * width + j] < 0)
			{
				image[i * width + j] = 0;
			}
			if (image[i * width + j] > 255)
			{
				image[i * width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;
}
int dis(int a, int x)
{
	return pow((a - x),2);
}
int assignCluster(int point, int* cluster_points, int K)
{
	int min = INT_MAX;
	int min_ind;
	for (int j = 0; j < K; j++)
	{
		int dist = dis(point, cluster_points[j]);
		if (dist < min)
		{
			min = dist;
			min_ind = j;
		}
	}
	return min_ind;
}

int main()
{
	//###############3MPI############################3
	// Initial MPI and find process rank and number of processes.
	cout << "STARTING.......................................\n";
	MPI_Init(NULL, NULL);
	int ImageWidth = 4, ImageHeight = 4;
	int start_s, stop_s, TotalTime = 0;

	System::String^ imagePath;
	std::string img;
	img = "..//Data//Input//test1.png";

	imagePath = marshal_as<System::String^>(img);

	start_s = clock();

	int rank, nprocs;
	MPI_Status stat;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	int tag = 200;
	nCentroids = 2;
	int K = nCentroids;


	int*data = inputImage(&ImageWidth, &ImageHeight, imagePath);
	IMAGE_BYTES = ImageWidth * ImageHeight * sizeof(int);
	CLUSTER_BYTES = nCentroids * sizeof(int);
	Byte_size = ImageWidth * ImageHeight;
	bytes_color = (int*)(malloc(IMAGE_BYTES));
	ByteCentroid = (int*)(malloc(CLUSTER_BYTES));
    for(int i=0;i< Byte_size;i++)
	         bytes_color[i] = data[i];

	int* ImageData = new int[Byte_size];
	int points_per_node = Byte_size / (nprocs - 1);
	//mpiexec "HPC_ProjectTemplate.exe"
	if (rank == 0) {

		int* assignment = new int[Byte_size];
		cout << "Random Centroids of  " << K;
		int* cluster_points = new int[K];
		// to get random centroid for each k
		for (int i = 0; i < K; i++)
		{
			int random_pt = rand() % Byte_size;
			cluster_points[i] = bytes_color[random_pt];
			cout << "\n"<<i << " :   "<<cluster_points[i];
		}
		// separte data on #procssor
		for (int iter = 0; iter < 100; iter++)
		{
			//cout << "ieration " << iter << endl;
			int start_index = 0, end_index = -1;
			//cout << "send processor/n";
			for (int j = 0; j < nprocs - 2; j++)
			{
				start_index = j * points_per_node;
				end_index = (j + 1) * points_per_node - 1;
				int num_indices = end_index - start_index + 1;
				MPI_Send(&start_index, 1, MPI_INT, j + 1, tag, MPI_COMM_WORLD);
				MPI_Send(&end_index, 1, MPI_INT, j + 1, tag, MPI_COMM_WORLD);
			}
			end_index += 1;
			MPI_Send(&end_index, 1, MPI_INT, nprocs - 1, tag, MPI_COMM_WORLD);
			int s = Byte_size;
			MPI_Send(&(s), 1, MPI_INT, nprocs - 1, tag, MPI_COMM_WORLD);

			//cout << "Sepate data/n";
			//MPI_Send(cluster_points, K, MPI_INT, 1, tag, MPI_COMM_WORLD);
			MPI_Bcast(cluster_points, K, MPI_INT, 0, MPI_COMM_WORLD);

			long * cluster_sum = new long[K];
			long * cluster_count = new long[K];
			memset(cluster_count, 0, sizeof(cluster_count));
			memset(cluster_sum, 0, sizeof(cluster_sum));

			//cout << "take data from processor/n";
			for (int proc = 1; proc < nprocs; proc++)
			{
				for (int i = 0; i < K; i++)
				{
					long sum, count;
					MPI_Recv(&sum, 1, MPI_LONG, proc, tag, MPI_COMM_WORLD, &stat);
					MPI_Recv(&count, 1, MPI_LONG, proc, tag, MPI_COMM_WORLD, &stat);
					//cout <<sum << "  , " << count << "\n";
					cluster_sum[i] += sum;
					cluster_count[i] += count;
					//cout << cluster_sum[i] << "  , " << cluster_count[i] << "\n";
				}
			}
			//cout << "clusters_points\n";
			for (int i = 0; i < K; i++)
			{
			//	cout << cluster_sum[i] << "  , " << cluster_count[i]<<"\n";
				cluster_points[i] =abs( cluster_sum[i] / cluster_count[i]);
			
			}

			if (iter == 99)
			{
				//cout << "Last step/n";
				int sti, endi;
				for (int i = 1; i < nprocs; i++) {
					//cout << "recev....";
					MPI_Recv(&sti, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &stat);
					MPI_Recv(&endi, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &stat);
					int* assignment_k = new int[endi - sti + 1];
					MPI_Recv(assignment_k, endi - sti + 1, MPI_INT, i, tag, MPI_COMM_WORLD, &stat);
					for (int j = sti; j <= endi; j++)
					{
						assignment[j] = assignment_k[j - sti];
					}
				}
				for (int i = 0; i < Byte_size; i++)
				{
					ImageData[i] = cluster_points[assignment[i]];
				}
			}
		}
		cout << "\nAnswer:";
		for (int i = 0; i < K; i++)
		{
			cout << "\n" << i << " :   " << cluster_points[i];
		}
		cout<<"\n";
		createImage(ImageData, ImageWidth, ImageHeight, 3);
		/*
		System::Drawing::Bitmap MyNewImage(ImageWidth, ImageHeight);
		
		for (int i = 0; i < MyNewImage.Height; i++)
		{
			for (int j = 0; j < MyNewImage.Width; j++)
			{
				System::Drawing::Color c = System::Drawing::Color::FromArgb(ImageData[i * ImageWidth + j], ImageData[i * ImageWidth + j], ImageData[i * ImageWidth + j]);
				MyNewImage.SetPixel(j, i, c);
			}
		}
		MyNewImage.Save("..//Data//Output//outputRes" + 3 + ".png");
		cout << "result Image Saved " << 1 << endl;*/
	}
	else {
		for (int iter = 0; iter < 100; iter++)
		{

			int start, end;
			MPI_Recv(&start, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &stat);
			MPI_Recv(&end, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &stat);

			int* cluster_points = new int[K];

			//MPI_Recv(cluster_points, K, MPI_INT, 0, tag, MPI_COMM_WORLD,&stat);
			MPI_Bcast(cluster_points, K, MPI_INT, 0, MPI_COMM_WORLD);

			//cout<<"Here already\n";
			//for(int i=0;i<K;i++)
			//	cout<<i<<" : "<<cluster_points[i]<<"  , ";
			//cout<<"\n";

			map<int, vector<int> > cluster;
			int* assignment = new int[end - start + 1];

			for (int i = start; i <= end; i++)
			{
				//cout << "data :" << bytes_color[i] << " , ";
				int c = assignCluster(bytes_color[i], cluster_points, K);
				// if exist in cluster
				if (cluster.find(c) == cluster.end())
				{
					vector<int> v;
					v.push_back(i);
					cluster.insert(make_pair(c, v));
				}
				else
					cluster[c].push_back(i);
				assignment[i - start] = c;
				//cout << "During " << c<<"\n";
			}
			for (int i = 0; i < K; i++)
			{   
				     long sum = 0, count = 0;
					for (int j = 0; j < cluster[i].size(); j++)
					{
						//cout << bytes_color[cluster[i][j]] <<"   ,    ";
						sum += bytes_color[cluster[i][j]];

					}
					count = cluster[i].size();
				
				MPI_Send(&sum, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
				MPI_Send(&count, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
			}
			if (iter == 99)
			{
				MPI_Send(&start, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
				MPI_Send(&end, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
				MPI_Send(assignment, end - start + 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
				
			}
		}
	}

	stop_s = clock();
	TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
	cout << "time: " << TotalTime << endl;
	//free(bytes_color);

	MPI_Finalize();
	return 0;

}



