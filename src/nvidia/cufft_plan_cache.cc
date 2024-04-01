/* Copyright (c) 2024 NVIDIA Corporation. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. */

#include <iostream>
#include <unordered_map>
typedef std::tuple<int, int, int, int, int> fft_param;

class PlanHashFunction {
public:
	size_t operator() (const fft_param& p) const {
		return std::get<0>(p) * std::get<4>(p) * (std::get<1>(p) + std::get<2>(p) + std::get<3>(p));
	}
};


std::unordered_map<fft_param, int, PlanHashFunction> plan_table;

extern "C" {
void find_fft_plan(int inc, int jump, int n, int m, int isign, int *plan) {
  fft_param param(inc, jump, n, m, isign);
  auto it = plan_table.find(param);
  if (it == plan_table.end()) {
    *plan = -1;
  } else {
    *plan = it->second;
  }
}

void cache_fft_plan(int inc, int jump, int n, int m, int isign, int plan) {
  fft_param param(inc, jump, n, m, isign);
  plan_table[param] = plan;
}

void fft_plan_size(int *size) { *size = plan_table.size(); }

void get_plan_list(int *plan_list) {
  int index = 0;
  for (auto it = plan_table.begin(); it != plan_table.end(); ++it) {
    *(plan_list + index) = it->second;
    ++index;
  }
}
}

