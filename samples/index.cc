#include <efanna.hpp>
#include <iostream>
#include <fstream>
#include <chrono>

using namespace efanna;
using namespace std;
void load_data(char* filename, float*& data, size_t& num,int& dim){// load data with sift10K pattern
  ifstream in(filename, ios::binary);
  if(!in.is_open()){cout<<"open file error"<<endl;exit(-1);}
  in.read((char*)&dim,4);
  cout<<"data dimension: "<<dim<<endl;
  in.seekg(0,ios::end);
  ios::pos_type ss = in.tellg();
  size_t fsize = (size_t)ss;
  num = fsize / (dim+1) / 4;
  data = new float[num*dim];

  in.seekg(0,ios::beg);
  for(size_t i = 0; i < num; i++){
    in.seekg(4,ios::cur);
    in.read((char*)(data+i*dim),dim*4);
  }
  in.close();
}



int main(int argc, char** argv){
  if(argc!=4){cout<< argv[0] << " data_file index_file tableNum)" <<endl; exit(-1);}

  float* data_load = NULL;
  size_t points_num;
  int dim;
  load_data(argv[1], data_load, points_num,dim);
  Matrix<float> dataset(points_num,dim,data_load);

  int table = atoi(argv[3]);
  if (table < 0 || table*32 > 100000){cout<<"tableNum error!";exit(-1);}


  FIndex<float> index(dataset, new L2DistanceAVX<float>(), efanna::HAMMINGIndexParams(table));

  auto s = std::chrono::high_resolution_clock::now();

  index.buildIndex(data_load, points_num, dim);

  auto e = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = e-s;
  std::cout << "indexing time: " << diff.count() << "\n";

  index.saveIndex(argv[2]);

  return 0;
}
