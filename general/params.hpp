#ifndef EFANNA_PARAMS_H_
#define EFANNA_PARAMS_H_
#include <map>
#include <string>

namespace efanna {
  enum init_algorithm {
      HAMMING,
  };

  union ValueType{
    int int_val;
    float float_val;
    char* str_pt;
  };


  typedef std::map<std::string, ValueType> ExtraParamsMap;
  struct IndexParams{
    init_algorithm init_index_type;
    unsigned index_tables=32;
    unsigned nGroup=1000;
    unsigned nIter=100;
    ExtraParamsMap extra_params;
  };

  struct SearchParams{
    unsigned search_init_num;
    unsigned search_groups;
  };

}
#endif
