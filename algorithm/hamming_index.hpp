// Copyright (C) 2018 Deng Cai <dengcai@gmail.com>. All Rights Reserved.


#ifndef EFANNA_HAMMING_INDEX_H_
#define EFANNA_HAMMING_INDEX_H_

#include "algorithm/base_index.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <cblas.h>



namespace efanna{


inline unsigned parallel_popcnt32(unsigned x) {
	x = (x & 0x55555555) + ((x >> 1 ) & 0x55555555);
	x = (x & 0x33333333) + ((x >> 2 ) & 0x33333333);
	x = (x & 0x0f0f0f0f) + ((x >> 4 ) & 0x0f0f0f0f);
	x = (x & 0x00ff00ff) + ((x >> 8 ) & 0x00ff00ff);
	x = (x & 0x0000ffff) + ((x >> 16) & 0x0000ffff);
	return x;
}

struct HAMMINGIndexParams : public IndexParams
{
	HAMMINGIndexParams(int nCluster, int iter, int table)
	{
		init_index_type = HAMMING;
		index_tables = table;
		nGroup=nCluster;
		nIter=iter;
	}
};

template <typename DataType>
class HAMMINGIndex : public InitIndex<DataType>
{
public:

	typedef InitIndex<DataType> BaseClass;

	typedef std::vector<unsigned> Code;
	typedef std::vector<Code> Codes;
	typedef std::vector<Codes> BucketCodes;


	HAMMINGIndex(const Matrix<DataType>& dataset, const Distance* d, const IndexParams& params = HAMMINGIndexParams(0)) :
		BaseClass(dataset,d,params)
	{}

	~HAMMINGIndex(){
		delete projection_matrix;
	}

	void loadIndex(char* filename){
		std::ifstream in(filename, std::ios::binary);
		if(!in.is_open()){std::cout<<"open file error"<<std::endl;exit(-1);}

		unsigned dim;
		unsigned codelen;

		in.read((char*)&dim,sizeof(int));
		in.read((char*)&codelen,sizeof(int));

		projection_matrix = new float[codelen * dim];
    	for (unsigned i = 0; i < codelen * dim; i++) {
    		in.read((char*)&projection_matrix[i], sizeof(int));
    	}

    	unsigned points_num;
    	unsigned nTables;
    	in.read((char*)&points_num, sizeof(int));
    	in.read((char*)&nTables, sizeof(int));

    	for(unsigned j = 0; j < nTables; j++){
    		Code b;
    		for(unsigned i = 0; i < points_num; i++){
    			unsigned int codetmp;
    			in.read((char*)&codetmp,sizeof(int));
    			b.push_back(codetmp);
    		}
    		BaseCode.push_back(b);
    	}

    	unsigned nCluster;
		in.read((char*)&nCluster,sizeof(int));
		centers.resize(nCluster);
		for(unsigned i = 0; i < nCluster; i++){
			for(unsigned j = 0; j < dim; j++){
				float tmp;
				in.read((char*)&tmp,sizeof(float));
				centers[i].push_back(tmp);
			}
		}

		for(unsigned i = 0; i < points_num; i++){
			unsigned indextmp;
			in.read((char*)&indextmp,sizeof(int));
			labelIndex.push_back(indextmp);
		}

		in.close();

		hashTable.resize(nCluster);
		for(unsigned i = 0; i < nCluster; i++){
			Codes tmp(nTables);
			BaseBucket.push_back(tmp);
		}

		for(unsigned i = 0; i < points_num; i++){
			hashTable[labelIndex[i]].push_back(i);
			for(unsigned j=0; j<nTables; j++){
				BaseBucket[labelIndex[i]][j].push_back(BaseCode[j][i]);
			}
		}



		std::cout << "Index loaded!" << std::endl;


    }

    void saveIndex(char* filename){
    	std::ofstream out(filename,std::ios::binary);
		if(!out.is_open()){std::cout<<"open file error"<<std::endl;exit(-1);}

		unsigned dim = features_.get_cols();
		unsigned codelen = params_.index_tables * 32;
    	out.write((char*)&dim, sizeof(int));
    	out.write((char*)&codelen, sizeof(int));

    	for (unsigned i = 0; i < codelen * dim; i++) {
    		out.write((char*)&projection_matrix[i], sizeof(int));
    	}

    	unsigned points_num = features_.get_rows();
    	out.write((char*)&points_num, sizeof(int));
    	out.write((char*)&params_.index_tables, sizeof(int));


		for (unsigned j = 0; j < params_.index_tables; j++) {
			for (unsigned i = 0; i < points_num; i++) {
				out.write((char*)&BaseCode[j][i], sizeof(int));
			}
		}

		unsigned nCluster = centers.size();
		out.write((char*)&nCluster, sizeof(int));
		for(unsigned i=0;i<nCluster;i++){
			for (unsigned j = 0; j < dim; j++) {
				float tmp=centers[i][j];
				out.write((char*)&tmp, sizeof(float));
			}
		}

		for (unsigned i = 0; i < points_num; i++) {
			out.write((char*)&labelIndex[i], sizeof(int));
		}

    	out.close();
    }


	void generate_random_projection_matrix(int dim, int codelen, float * projection_matrix) {
	  std::default_random_engine generator;
	  std::normal_distribution<float> distribution(0.0, 1.0);
	  for (int i = 0; i < codelen * dim; i++) {
	    projection_matrix[i] = distribution(generator);
	  }
	}

	void random_projection(float *data, size_t points_num, int dim, int codelen,Codes &output) {
	  const int numTable = params_.index_tables;
	  float *projection_result = new float[codelen * points_num];
	  cblas_sgemm(CblasColMajor, CblasTrans, CblasTrans, points_num, codelen, dim,
	              1, data, dim, projection_matrix, codelen, 0, projection_result, points_num);

	  output.resize(numTable);
	  for (int i = 0 ; i < numTable; i++) {
		  output[i].resize(points_num);
	    std::fill(output[i].begin(), output[i].end(), 0);
	  }

	  for (int i = 0 ; i < codelen; i++) {
	    const int table_id = i / 32;
	    const unsigned bit = 1u << (i % 32);
	    for (size_t j = 0; j < points_num; j++) {
	      if (projection_result[i * points_num + j] > 0)
	    	  output[table_id][j] |= bit;
	      }
	  }

	  delete []projection_result;
	}

	void kmeansPartition_fast(){
		unsigned nCluster = params_.nGroup;

		unsigned nSmp = features_.get_rows();
		unsigned mFea = features_.get_cols();

		DistanceFastL2 *dist_ = (DistanceFastL2*)distance_;

		std::vector<unsigned> labelIndexNew(nSmp);
		std::vector<float> centerDistance(nSmp);
		std::vector<unsigned> pointNum(nCluster);
		std::vector<float> centerNorm(nCluster);

		labelIndex.resize(nSmp,0);
		centers.resize(nCluster);
		for(unsigned i = 0; i<nCluster; i++){
			centers[i].resize(mFea,0);
		}

		unsigned avg = nSmp / nCluster;
		for(unsigned i = 0; i<nCluster-1; i++){
			for(unsigned j = i*avg; j<i*avg+avg; j++){
				for (unsigned k=0; k<mFea; k++){
					centers[i][k] += features_.get_row(j)[k];
				}
			}
			for (unsigned k=0; k<mFea; k++){
				centers[i][k] = centers[i][k]/avg;
			}
		}
		for(unsigned j = (nCluster-1)*avg; j<nSmp; j++){
			for (unsigned k=0; k<mFea; k++){
				centers[nCluster-1][k] += features_.get_row(j)[k];
			}
		}
		for (unsigned k=0; k<mFea; k++){
			centers[nCluster-1][k] = centers[nCluster-1][k]/(nSmp-(nCluster-1)*avg+1);
		}

		for(unsigned i = 0; i<nCluster; i++){
			centerNorm[i]= dist_->norm(centers[i].data(),mFea);
		}


		unsigned iter=0;
		bool converged=false;
		while(iter<params_.nIter && !converged){

			#pragma omp parallel for
			for(unsigned i = 0; i<nSmp; i++){
				std::vector<float> dist(nCluster);
				for(unsigned k=0;k<nCluster;k++){
					dist[k] = dist_->compare(features_.get_row(i), centers[k].data(), centerNorm[k], mFea);
				}
				labelIndexNew[i]=std::min_element(dist.begin(),dist.end())-dist.begin();
				centerDistance[i] = dist[labelIndexNew[i]];
			}

			if(iter%10==0){
				std::cout <<"Iteration " <<iter<<" : "<< std::accumulate(centerDistance.begin(),centerDistance.end(),0.0) <<std::endl;
			}


			unsigned j=0;
			for(;j<nSmp;j++){
				if(labelIndex[j] != labelIndexNew[j]) break;
			}
			if(j==nSmp){
				converged = true;
			}else{
				iter++;

				for(unsigned i = 0; i<nCluster; i++){
					std::fill(centers[i].begin(), centers[i].end(), 0);
				}
				std::fill(pointNum.begin(), pointNum.end(), 0);

				for(unsigned i = 0; i<nSmp; i++){
					for (unsigned k=0; k<mFea; k++){
						centers[labelIndexNew[i]][k] += features_.get_row(i)[k];
					}
					pointNum[labelIndexNew[i]] ++;
				}

				// initialize index
				std::vector<unsigned> clusterIndex(nCluster);
				for (size_t i = 0; i < clusterIndex.size(); i++) {
					clusterIndex[i] = i;
				}
				// compare function
				auto compare_func = [&pointNum](const size_t a, const size_t b)->bool
						{
					return pointNum[a] < pointNum[b];
						};
				// sort the index by its value
				std::sort(clusterIndex.begin(), clusterIndex.end(),  compare_func);

				if(pointNum[clusterIndex[0]] == 0){
					std::cout <<"Empty cluster!"<<std::endl;
				}

				unsigned j=0;
				while(pointNum[clusterIndex[j]] == 0){
					unsigned mv=pointNum[clusterIndex[nCluster-1-j]]/2;
					unsigned mved=0;
					for(unsigned i = 0; i<nSmp; i++){
						if(labelIndexNew[i] == clusterIndex[nCluster-1-j] && mved<mv){
							for (unsigned k=0; k<mFea; k++){
								centers[clusterIndex[j]][k] += features_.get_row(i)[k];
								centers[clusterIndex[nCluster-1-j]][k] -= features_.get_row(i)[k];
							}
							pointNum[clusterIndex[j]]++;
							pointNum[clusterIndex[nCluster-1-j]]--;
							mved++;
						}
					}
					j++;
				}

				for(unsigned i = 0; i<nCluster; i++){
					for (unsigned k=0; k<mFea; k++){
						centers[i][k] = centers[i][k]/pointNum[i];
					}
				}
				for(unsigned i = 0; i<nCluster; i++){
					centerNorm[i]= dist_->norm(centers[i].data(),mFea);
				}
				labelIndex.swap(labelIndexNew);
			}
		}

		std::cout<<"Converged: "<<converged<<std::endl;



	}

	void kmeansPartition(){
		unsigned nCluster = params_.nGroup;

		unsigned nSmp = features_.get_rows();
		unsigned mFea = features_.get_cols();

		std::vector<unsigned> labelIndexNew(nSmp);
		std::vector<float> centerDistance(nSmp);
		std::vector<unsigned> pointNum(nCluster);

		labelIndex.resize(nSmp,0);
		centers.resize(nCluster);
		for(unsigned i = 0; i<nCluster; i++){
			centers[i].resize(mFea,0);
		}

		unsigned avg = nSmp / nCluster;
		for(unsigned i = 0; i<nCluster-1; i++){
			for(unsigned j = i*avg; j<i*avg+avg; j++){
				for (unsigned k=0; k<mFea; k++){
					centers[i][k] += features_.get_row(j)[k];
				}
			}
			for (unsigned k=0; k<mFea; k++){
				centers[i][k] = centers[i][k]/avg;
			}
		}
		for(unsigned j = (nCluster-1)*avg; j<nSmp; j++){
			for (unsigned k=0; k<mFea; k++){
				centers[nCluster-1][k] += features_.get_row(j)[k];
			}
		}
		for (unsigned k=0; k<mFea; k++){
			centers[nCluster-1][k] = centers[nCluster-1][k]/(nSmp-(nCluster-1)*avg+1);
		}

		unsigned iter=0;
		bool converged=false;
		while(iter<params_.nIter && !converged){

			#pragma omp parallel for
			for(unsigned i = 0; i<nSmp; i++){
				std::vector<float> dist(nCluster);
				for(unsigned k=0;k<nCluster;k++){
					dist[k] = distance_->compare(features_.get_row(i), centers[k].data(), mFea);
				}
				labelIndexNew[i]=std::min_element(dist.begin(),dist.end())-dist.begin();
				centerDistance[i] = dist[labelIndexNew[i]];
			}

			if(iter%10==0){
				std::cout <<"Iteration " <<iter<<" : "<< std::accumulate(centerDistance.begin(),centerDistance.end(),0.0) <<std::endl;
			}


			unsigned j=0;
			for(;j<nSmp;j++){
				if(labelIndex[j] != labelIndexNew[j]) break;
			}
			if(j==nSmp){
				converged = true;
			}else{
				iter++;

				for(unsigned i = 0; i<nCluster; i++){
					std::fill(centers[i].begin(), centers[i].end(), 0);
				}
				std::fill(pointNum.begin(), pointNum.end(), 0);

				for(unsigned i = 0; i<nSmp; i++){
					for (unsigned k=0; k<mFea; k++){
						centers[labelIndexNew[i]][k] += features_.get_row(i)[k];
					}
					pointNum[labelIndexNew[i]] ++;
				}

				// initialize index
				std::vector<unsigned> clusterIndex(nCluster);
				for (size_t i = 0; i < clusterIndex.size(); i++) {
					clusterIndex[i] = i;
				}
				// compare function
				auto compare_func = [&pointNum](const size_t a, const size_t b)->bool
						{
					return pointNum[a] < pointNum[b];
						};
				// sort the index by its value
				std::sort(clusterIndex.begin(), clusterIndex.end(),  compare_func);

				if(pointNum[clusterIndex[0]] == 0){
					std::cout <<"Empty cluster!"<<std::endl;
				}

				unsigned j=0;
				while(pointNum[clusterIndex[j]] == 0){
					unsigned mv=pointNum[clusterIndex[nCluster-1-j]]/2;
					unsigned mved=0;
					for(unsigned i = 0; i<nSmp; i++){
						if(labelIndexNew[i] == clusterIndex[nCluster-1-j] && mved<mv){
							for (unsigned k=0; k<mFea; k++){
								centers[clusterIndex[j]][k] += features_.get_row(i)[k];
								centers[clusterIndex[nCluster-1-j]][k] -= features_.get_row(i)[k];
							}
							pointNum[clusterIndex[j]]++;
							pointNum[clusterIndex[nCluster-1-j]]--;
							mved++;
						}
					}
					j++;
				}

				for(unsigned i = 0; i<nCluster; i++){
					for (unsigned k=0; k<mFea; k++){
						centers[i][k] = centers[i][k]/pointNum[i];
					}
				}
				labelIndex.swap(labelIndexNew);
			}
		}

		std::cout<<"Converged: "<<converged<<std::endl;



	}

	/*void kmeansPartitionImpl(const Matrix<DataType>& fea, std::vector<std::vector<float>>& cent, unsigned nIter){

		unsigned nSmp = features.get_rows();
		unsigned mFea = features.get_cols();

		std::vector<unsigned> labelIndexNew(nSmp);
		std::vector<float> centerDistance(nSmp);
		std::vector<unsigned> pointNum(nCluster);

		labelIndex.resize(nSmp,0);
		centers.resize(nCluster);
		for(unsigned i = 0; i<nCluster; i++){
			centers[i].resize(mFea,0);
		}

		unsigned avg = nSmp / nCluster;
		for(unsigned i = 0; i<nCluster-1; i++){
			for(unsigned j = i*avg; j<i*avg+avg; j++){
				for (unsigned k=0; k<mFea; k++){
					centers[i][k] += features_.get_row(j)[k];
				}
			}
			for (unsigned k=0; k<mFea; k++){
				centers[i][k] = centers[i][k]/avg;
			}
		}
		for(unsigned j = (nCluster-1)*avg; j<nSmp; j++){
			for (unsigned k=0; k<mFea; k++){
				centers[nCluster-1][k] += features_.get_row(j)[k];
			}
		}
		for (unsigned k=0; k<mFea; k++){
			centers[nCluster-1][k] = centers[nCluster-1][k]/(nSmp-(nCluster-1)*avg+1);
		}

		unsigned iter=0;
		bool converged=false;
		while(iter<params_.nIter && !converged){

			#pragma omp parallel for
			for(unsigned i = 0; i<nSmp; i++){
				std::vector<float> dist(nCluster);
				for(unsigned k=0;k<nCluster;k++){
					dist[k] = distance_->compare(features_.get_row(i), centers[k].data(), mFea);
				}
				labelIndexNew[i]=std::min_element(dist.begin(),dist.end())-dist.begin();
				centerDistance[i] = dist[labelIndexNew[i]];
			}

			if(iter%10==0){
				std::cout <<"Iteration " <<iter<<" : "<< std::accumulate(centerDistance.begin(),centerDistance.end(),0.0) <<std::endl;
			}


			unsigned j=0;
			for(;j<nSmp;j++){
				if(labelIndex[j] != labelIndexNew[j]) break;
			}
			if(j==nSmp){
				converged = true;
			}else{
				iter++;

				for(unsigned i = 0; i<nCluster; i++){
					std::fill(centers[i].begin(), centers[i].end(), 0);
				}
				std::fill(pointNum.begin(), pointNum.end(), 0);

				for(unsigned i = 0; i<nSmp; i++){
					for (unsigned k=0; k<mFea; k++){
						centers[labelIndexNew[i]][k] += features_.get_row(i)[k];
					}
					pointNum[labelIndexNew[i]] ++;
				}

				// initialize index
				std::vector<unsigned> clusterIndex(nCluster);
				for (size_t i = 0; i < clusterIndex.size(); i++) {
					clusterIndex[i] = i;
				}
				// compare function
				auto compare_func = [&pointNum](const size_t a, const size_t b)->bool
						{
					return pointNum[a] < pointNum[b];
						};
				// sort the index by its value
				std::sort(clusterIndex.begin(), clusterIndex.end(),  compare_func);

				if(pointNum[clusterIndex[0]] == 0){
					std::cout <<"Empty cluster!"<<std::endl;
				}

				unsigned j=0;
				while(pointNum[clusterIndex[j]] == 0){
					unsigned mv=pointNum[clusterIndex[nCluster-1-j]]/2;
					unsigned mved=0;
					for(unsigned i = 0; i<nSmp; i++){
						if(labelIndexNew[i] == clusterIndex[nCluster-1-j] && mved<mv){
							for (unsigned k=0; k<mFea; k++){
								centers[clusterIndex[j]][k] += features_.get_row(i)[k];
								centers[clusterIndex[nCluster-1-j]][k] -= features_.get_row(i)[k];
							}
							pointNum[clusterIndex[j]]++;
							pointNum[clusterIndex[nCluster-1-j]]--;
							mved++;
						}
					}
					j++;
				}

				for(unsigned i = 0; i<nCluster; i++){
					for (unsigned k=0; k<mFea; k++){
						centers[i][k] = centers[i][k]/pointNum[i];
					}
				}
				labelIndex.swap(labelIndexNew);
			}
		}

		std::cout<<"Converged: "<<converged<<std::endl;

	}*/

	void buildIndexImpl(float* data, size_t points_num, int dim){
		int codelen = params_.index_tables *32;
		projection_matrix = new float[codelen * dim];
		generate_random_projection_matrix(dim, codelen, projection_matrix);
		random_projection(data, points_num, dim, codelen, BaseCode);
		kmeansPartition_fast();
	}

	void getNeighbors(size_t K, const Matrix<DataType>& query, float* query_){

		unsigned nCluster = centers.size();
		unsigned nTable = BaseCode.size();

		unsigned nGroup = SP.search_groups;

		if(nGroup > nCluster){
			nGroup = nCluster;
		}

		random_projection(query_, query.get_rows(), query.get_cols(), nTable*32, QueryCode);

		std::vector<std::vector<int> >hammingDistance(nGroup);
		std::vector<size_t> hammingDistanceIndex;;
		std::vector<std::pair<int, unsigned> >hammingDistanceSort;

		nn_results.clear();
		for (size_t cur = 0; cur < query.get_rows(); cur++) {


			std::vector<std::pair<float, unsigned> > CenterDistance;
			for(unsigned j = 0; j < nCluster; j++){
				CenterDistance.push_back(std::make_pair(distance_->compare(query.get_row(cur), centers[j].data(), query.get_cols()), j));
			}
			std::partial_sort(CenterDistance.begin(), CenterDistance.begin() + nGroup, CenterDistance.end());

			size_t pNum = 0;
			std::vector<unsigned> compareIndex;

			for(unsigned j = 0; j < nGroup; j++){
				compareIndex.push_back(CenterDistance[j].second);
				pNum += hashTable[CenterDistance[j].second].size();

			}
			if ((unsigned)SP.search_init_num > pNum){
				SP.search_init_num = pNum;
			}


			// resize hammingDistance
			for (size_t i = 0; i < compareIndex.size(); i++) {
				hammingDistance[i].resize(hashTable[compareIndex[i]].size());
				std::fill(hammingDistance[i].begin(), hammingDistance[i].end(), 0);
			}

			unsigned IndexNum = compareIndex.size();
			for (unsigned j = 0; j < nTable; j++) {
				for (unsigned idx = 0; idx < IndexNum; idx++) {
					unsigned tNum = BaseBucket[compareIndex[idx]][j].size();
					Code &code = BaseBucket[compareIndex[idx]][j];
					std::vector<int> &dst = hammingDistance[idx];
					for (unsigned i = 0; i < tNum; i++) {
						dst[i] += parallel_popcnt32(code[i] ^ QueryCode[j][cur]);
					}
				}
			}


			hammingDistanceSort.resize(0);
			for (unsigned idx = 0; idx < IndexNum; idx++) {
				size_t sort_num = hammingDistance[idx].size();
				for (size_t i = 0; i < sort_num; i++) {
					hammingDistanceSort.push_back(std::pair<int, unsigned>(hammingDistance[idx][i], hashTable[compareIndex[idx]][i]));
				}
			}



			// sort the index by its value
			size_t sort_num = std::min((size_t)SP.search_init_num, hammingDistanceSort.size());
			std::partial_sort(hammingDistanceSort.begin(), hammingDistanceSort.begin() + sort_num, hammingDistanceSort.end());

			std::vector<int> res;
			if ((unsigned) SP.search_init_num <= K){
				for (unsigned int j = 0; j < (unsigned) sort_num; j++) res.push_back(hammingDistanceSort[j].second);
			}else{
				std::vector<std::pair<float, unsigned> > Distance;
				for (size_t i = 0; i < (unsigned) sort_num; i++) {
					Distance.push_back(std::make_pair(distance_->compare(query.get_row(cur), features_.get_row(hammingDistanceSort[i].second), features_.get_cols()), hammingDistanceSort[i].second));
				}
				std::partial_sort(Distance.begin(), Distance.begin() + K, Distance.end());

				for (unsigned int j = 0; j < K; j++) res.push_back(Distance[j].second);
			}
			nn_results.push_back(res);
		}
	}


protected:
	int codelength;
	float *projection_matrix = NULL;

	std::vector<std::vector<float>> centers;
	std::vector<unsigned> labelIndex;

	std::vector<std::vector<unsigned>> hashTable;
	BucketCodes BaseBucket;


	Codes BaseCode;
	Codes QueryCode;

	USING_BASECLASS_SYMBOLS

};
}
#endif
