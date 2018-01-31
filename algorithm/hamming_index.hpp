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
	HAMMINGIndexParams(int table)
	{
		init_index_type = HAMMING;
		index_tables = table;
	}
};

template <typename DataType>
class HAMMINGIndex : public InitIndex<DataType>
{
public:

	typedef InitIndex<DataType> BaseClass;

	typedef std::vector<unsigned> Code;
	typedef std::vector<Code> Codes;


	HAMMINGIndex(const Matrix<DataType>& dataset, const Distance<DataType>* d, const IndexParams& params = HAMMINGIndexParams(0)) :
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

		in.close();

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


		for (size_t j = 0; j < params_.index_tables; j++) {
			for (size_t i = 0; i < points_num; i++) {
				out.write((char*)&BaseCode[j][i], sizeof(int));
			}
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

	void buildIndexImpl(float* data, size_t points_num, int dim){
		int codelen = params_.index_tables *32;
		projection_matrix = new float[codelen * dim];
		generate_random_projection_matrix(dim, codelen, projection_matrix);
		random_projection(data, points_num, dim, codelen, BaseCode);

	}

	void getNeighbors(size_t K, const Matrix<DataType>& query, float* query_){

		size_t pNum = BaseCode[0].size();
		size_t nTable = BaseCode.size();

		random_projection(query_, query.get_rows(), query.get_cols(), nTable*32, QueryCode);


		if (SP.search_tables > nTable){
			std::cout << "The index only contains "<< nTable <<" tables, so "<<nTable<<" used."<<std::endl;
		}else{
			nTable =SP.search_tables;
		}


		nn_results.clear();
		for (size_t cur = 0; cur < query.get_rows(); cur++) {
			std::vector<int>hammingDistance(pNum);
			std::fill(hammingDistance.begin(), hammingDistance.end(), 0);

			for (size_t j = 0; j < nTable; j++) {
				for (size_t i = 0; i < pNum; i++) {
					hammingDistance[i] += parallel_popcnt32(BaseCode[j][i] ^ QueryCode[j][cur]);
				}
			}

			// initialize index
			std::vector<size_t> hammingDistanceIndex(pNum);
			for (size_t i = 0; i < hammingDistanceIndex.size(); i++) {
				hammingDistanceIndex[i] = i;
			}
			// compare function
			auto compare_func = [&hammingDistance](const size_t a, const size_t b)->bool
					{
				return hammingDistance[a] < hammingDistance[b];
					};
			// sort the index by its value
			std::partial_sort(hammingDistanceIndex.begin(), hammingDistanceIndex.begin() + SP.search_init_num, hammingDistanceIndex.end(),  compare_func);

			std::vector<std::pair<float, size_t>> Distance;
			for (size_t i = 0; i < (unsigned) SP.search_init_num; i++) {
				Distance.push_back(std::make_pair(distance_->compare(query.get_row(cur), features_.get_row(hammingDistanceIndex[i]), features_.get_cols()), hammingDistanceIndex[i]));
			}
			std::partial_sort(Distance.begin(), Distance.begin() + K, Distance.end());

			std::vector<int> res;
			for (unsigned int j = 0; j < K; j++) res.push_back(Distance[j].second);
			nn_results.push_back(res);
		}


	}


protected:
	int codelength;
	float *projection_matrix = NULL;

	Codes BaseCode;
	Codes QueryCode;

	USING_BASECLASS_SYMBOLS

};
}
#endif
