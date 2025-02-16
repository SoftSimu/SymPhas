
#include "testvirtualmodel.h"

#include <filesystem>
#include <fstream>

#include "symphas.h"

void testvirtualmodel() {
  // These two files must exist and be in checkpoint format.
  const char* first_file = "tests\\data0";
  const char* second_file = "tests\\data02";

  const char* test_directory = "tests";
  for (iter_type i = 0; i < 10; ++i) {
    char name[64]{};
    symphas::io::copy_data_file_name(test_directory, "", i, i,
                                     DataFileType::CHECKPOINT_DATA, name);

    FILE* f;
    if ((f = fopen(name, "w")) == 0) {
      symphas::lib::make_directory_for_file(name);
      f = fopen(name, "w");
    }
    fprintf(f, "%d\n", i);
    fclose(f);
  }

  char dir_with_name[256]{};
  symphas::io::copy_data_file_name(
      test_directory, "", 0, 0, DataFileType::CHECKPOINT_DATA, dir_with_name);
  const char* name_it = std::strrchr(dir_with_name, '/') + 1;
  char dir[256]{};
  char* dir_ptr = &dir[0];
  symphas::lib::get_parent_directory(dir_with_name, dir_ptr);

  for (auto& p : std::filesystem::directory_iterator(dir)) {
    char name[256]{};
    std::strcpy(name, p.path().string().c_str());
    const char* it = std::strrchr(name, '/');
    printf("%s, %d\n", it + 1, std::strcmp(it + 1, name_it));
  }

  // symphas::problem_parameters_type pp(1);

  // symphas::init_data_type tdataf{ { Axis::NONE, { Inside::FILE,
  // symphas::init_data_read{ first_file, 0 } } } };
  // pp.set_initial_data(&tdataf);

  SaveParams pyma_save{100, 1000};
  ModelVirtualBasic<2, scalar_t> pyma("load_dir", pyma_save);
  Conf c(symphas::conf::config());
  c.directory_settings.set_directory("pyma");
  c.simulation_settings.save = pyma_save;
  symphas::find_solution(pyma, c);

  // symphas::init_data_type tdata{ { Axis::NONE, { Inside::FILE, { second_file,
  // 1000 } } } }; System<scalar_t, 2> sys(tdata, pp.get_interval_data()[0]);
  //
  // auto mypy_end = sys.as_field();
  // auto sim_end = pyma.template system<0>().as_field();
  // symphas::io::save_vec(symphas::lib::abs(mypy_end - sim_end));
}
