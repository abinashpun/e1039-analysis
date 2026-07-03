#ifndef _UTIL_MINE__H_
#define _UTIL_MINE__H_
#include <string>
#include <vector>

/// Collection of your useful functions.
/*
 * Commonly-useful functions should be imported into e1039-core/packages/UtilAna.
 */
namespace UtilMine {
  std::string GetDstFilePath(const int run_id, const int spill_id, const std::string dir_dst="");
  std::vector<std::string> GetListOfSpillDSTs(const int run, const std::string dir_dst="");
};

#endif // _UTIL_MINE__H_
