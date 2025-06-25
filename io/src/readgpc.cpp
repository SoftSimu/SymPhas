
/* ***************************************************************************
 * This file is part of the SymPhas library, a framework for implementing
 * solvers for phase-field problems with compile-time symbolic algebra.
 *
 * Copyright (c) 2018-2021 by Steven A. Silber and Mikko Karttunen
 *
 * SymPhas is free software, which can be redistributed or modified under
 * the terms of the GNU Lesser General Public License (LGPL) as published
 * by the Free Software Foundation; LGPL version 3, or later versions at
 * your choice.
 *
 * SymPhas is distributed with the faith that it will be helpful and
 * practical but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * ***************************************************************************
 */

#include "readgpc.h"

template <>
void symphas::io::gp::col::read_block(scalar_t* grid,
                                      symphas::io::block_info binfo, FILE* f) {
  auto helper = symphas::io::new_helper(binfo);

  for (iter_type k = 0; k < GP_HELPER_LENZ; k++) {
    for (iter_type j = 0; j < GP_HELPER_LENY; j++) {
      for (iter_type i = 0; i < GP_HELPER_LENX; i++) {
        for (len_type n = 0; n < binfo.dimension(); ++n) {
          if (fscanf(f, "%*f") != 0) {
            fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
            exit(ERR_CODE_FILE_READ);
          }
        }
        scalar_t value;
        if (fscanf(f, "%lf", &value) != 1) {
          fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
          exit(ERR_CODE_FILE_READ);
        }
        if (grid != nullptr && helper->in_bounds({i, j, k})) {
          iter_type ii = GP_HELPER_INDEX({i, j, k});
          grid[ii] = value;
        }
      }
    }
  }
  symphas::io::free_helper(helper);
}

template <>
void symphas::io::gp::col::read_block(complex_t* grid,
                                      symphas::io::block_info binfo, FILE* f) {
  auto helper = symphas::io::new_helper(binfo);

  for (iter_type k = 0; k < GP_HELPER_LENZ; k++) {
    for (iter_type j = 0; j < GP_HELPER_LENY; j++) {
      for (iter_type i = 0; i < GP_HELPER_LENX; i++) {
        for (len_type n = 0; n < binfo.dimension(); ++n) {
          if (fscanf(f, "%*f") != 0) {
            fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
            exit(ERR_CODE_FILE_READ);
          }
        }
        double re, im;
        if (fscanf(f, "%lf+i%lf", &re, &im) == 0) {
          fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
          exit(ERR_CODE_FILE_READ);
        }

        if (grid != nullptr && helper->in_bounds({i, j, k})) {
          iter_type ii = GP_HELPER_INDEX({i, j, k});
          grid[ii] = complex_t{re, im};
        }
      }
    }
  }
  symphas::io::free_helper(helper);
}

template <>
void symphas::io::gp::col::read_block(double_arr2* grid,
                                      symphas::io::block_info binfo, FILE* f) {
  auto helper = symphas::io::new_helper(binfo);

  for (iter_type k = 0; k < GP_HELPER_LENZ; k++) {
    for (iter_type j = 0; j < GP_HELPER_LENY; j++) {
      for (iter_type i = 0; i < GP_HELPER_LENX; i++) {
        for (len_type n = 0; n < binfo.dimension(); ++n) {
          if (fscanf(f, "%*f") != 0) {
            fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
            exit(ERR_CODE_FILE_READ);
          }
        }

        double a, b;
        if (fscanf(f, "%lf %lf", &a, &b) == 0) {
          fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
          exit(ERR_CODE_FILE_READ);
        }
        if (grid != nullptr && helper->in_bounds({i, j, k})) {
          iter_type ii = GP_HELPER_INDEX({i, j, k});
          grid[ii][0] = a;
          grid[ii][1] = b;
        }
      }
    }
  }
  symphas::io::free_helper(helper);
}

template <>
void symphas::io::gp::col::read_block(vector_t<3>* grid,
                                      symphas::io::block_info binfo, FILE* f) {
  auto helper = symphas::io::new_helper(binfo);

  for (iter_type k = 0; k < GP_HELPER_LENZ; k++) {
    for (iter_type j = 0; j < GP_HELPER_LENY; j++) {
      for (iter_type i = 0; i < GP_HELPER_LENX; i++) {
        for (len_type n = 0; n < binfo.dimension(); ++n) {
          if (fscanf(f, "%*f") != 0) {
            fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
            exit(ERR_CODE_FILE_READ);
          }
        }

        double m, dx, dy, dz;
        if (fscanf(f, "%lf %lf %lf %lf", &dx, &dy, &dz, &m) == 0) {
          fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
          exit(ERR_CODE_FILE_READ);
        }

        if (grid != nullptr && helper->in_bounds({i, j, k})) {
          iter_type ii = GP_HELPER_INDEX({i, j, k});
          grid[ii] = vector_t<3>{dx * m, dy * m, dz * m};
        }
      }
    }
  }
  symphas::io::free_helper(helper);
}

template <>
void symphas::io::gp::col::read_block(vector_t<2>* grid,
                                      symphas::io::block_info binfo, FILE* f) {
  auto helper = symphas::io::new_helper(binfo);

  for (iter_type k = 0; k < GP_HELPER_LENZ; k++) {
    for (iter_type j = 0; j < GP_HELPER_LENY; j++) {
      for (iter_type i = 0; i < GP_HELPER_LENX; i++) {
        for (len_type n = 0; n < binfo.dimension(); ++n) {
          if (fscanf(f, "%*f") != 0) {
            fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
            exit(ERR_CODE_FILE_READ);
          }
        }

        double m, dx, dy;
        if (fscanf(f, "%lf %lf %lf", &dx, &dy, &m) == 0) {
          fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
          exit(ERR_CODE_FILE_READ);
        }

        if (grid != nullptr && helper->in_bounds({i, j, k})) {
          iter_type ii = GP_HELPER_INDEX({i, j, k});
          grid[ii] = vector_t<2>{dx * m, dy * m};
        }
      }
    }
  }
  symphas::io::free_helper(helper);
}

template <>
void symphas::io::gp::col::read_block(vector_t<1>* grid,
                                      symphas::io::block_info binfo, FILE* f) {
  auto helper = symphas::io::new_helper(binfo);

  for (iter_type k = 0; k < GP_HELPER_LENZ; k++) {
    for (iter_type j = 0; j < GP_HELPER_LENY; j++) {
      for (iter_type i = 0; i < GP_HELPER_LENX; i++) {
        for (len_type n = 0; n < binfo.dimension(); ++n) {
          if (fscanf(f, "%*f") != 0) {
            fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
            exit(ERR_CODE_FILE_READ);
          }
        }
        double m;
        if (fscanf(f, "%lf ", &m) == 0) {
          fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
          exit(ERR_CODE_FILE_READ);
        }

        if (grid != nullptr && helper->in_bounds({i, j, k})) {
          iter_type ii = GP_HELPER_INDEX({i, j, k});
          grid[ii] = vector_t<1>{m};
        }
      }
    }
  }
  symphas::io::free_helper(helper);
}

template <>
void symphas::io::gp::col::read_block(scalar_ptr_t (&grid)[3],
                                      symphas::io::block_info binfo, FILE* f) {
  auto helper = symphas::io::new_helper(binfo);

  for (iter_type k = 0; k < GP_HELPER_LENZ; k++) {
    for (iter_type j = 0; j < GP_HELPER_LENY; j++) {
      for (iter_type i = 0; i < GP_HELPER_LENX; i++) {
        for (len_type n = 0; n < binfo.dimension(); ++n) {
          if (fscanf(f, "%*f") != 0) {
            fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
            exit(ERR_CODE_FILE_READ);
          }
        }

        double m, dx, dy, dz;
        if (fscanf(f, "%lf %lf %lf %lf", &dx, &dy, &dz, &m) == 0) {
          fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
          exit(ERR_CODE_FILE_READ);
        }

        if (grid != nullptr && helper->in_bounds({i, j, k})) {
          iter_type ii = GP_HELPER_INDEX({i, j, k});
          grid[0][ii] = dx * m;
          grid[1][ii] = dy * m;
          grid[2][ii] = dz * m;
        }
      }
    }
  }
  symphas::io::free_helper(helper);
}

template <>
void symphas::io::gp::col::read_block(scalar_ptr_t (&grid)[2],
                                      symphas::io::block_info binfo, FILE* f) {
  auto helper = symphas::io::new_helper(binfo);

  for (iter_type k = 0; k < GP_HELPER_LENZ; k++) {
    for (iter_type j = 0; j < GP_HELPER_LENY; j++) {
      for (iter_type i = 0; i < GP_HELPER_LENX; i++) {
        for (len_type n = 0; n < binfo.dimension(); ++n) {
          if (fscanf(f, "%*f") != 0) {
            fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
            exit(ERR_CODE_FILE_READ);
          }
        }

        double m, dx, dy;
        if (fscanf(f, "%lf %lf %lf", &dx, &dy, &m) == 0) {
          fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
          exit(ERR_CODE_FILE_READ);
        }

        if (grid != nullptr && helper->in_bounds({i, j, k})) {
          iter_type ii = GP_HELPER_INDEX({i, j, k});
          grid[0][ii] = dx * m;
          grid[1][ii] = dy * m;
        }
      }
    }
  }
  symphas::io::free_helper(helper);
}

template <>
void symphas::io::gp::col::read_block(scalar_ptr_t (&grid)[1],
                                      symphas::io::block_info binfo, FILE* f) {
  auto helper = symphas::io::new_helper(binfo);

  for (iter_type k = 0; k < GP_HELPER_LENZ; k++) {
    for (iter_type j = 0; j < GP_HELPER_LENY; j++) {
      for (iter_type i = 0; i < GP_HELPER_LENX; i++) {
        for (len_type n = 0; n < binfo.dimension(); ++n) {
          if (fscanf(f, "%*f") != 0) {
            fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
            exit(ERR_CODE_FILE_READ);
          }
        }
        double m;
        if (fscanf(f, "%lf ", &m) == 0) {
          fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_READ, "scalar_t");
          exit(ERR_CODE_FILE_READ);
        }

        if (grid != nullptr && helper->in_bounds({i, j, k})) {
          iter_type ii = GP_HELPER_INDEX({i, j, k});
          grid[0][ii] = m;
        }
      }
    }
  }
  symphas::io::free_helper(helper);
}
