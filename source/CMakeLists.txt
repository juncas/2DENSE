# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

# 项目信息
project (2DENSE)
# 启动FORTRAN语言
enable_language(Fortran)

# 查找是否安装
find_package(OpenMP)
set(CMAKE_Fortran_COMPILER "gfortran")

set(CMAKE_Fortran_FLAGS "-ffree-line-length-200 -O2 ${OpenMP_Fortran_FLAGS}")

# 将所有FORTRAN文件插入
file(GLOB_RECURSE src_files ./*.f90)

# 生成目标文件
add_executable("2DENSE.out" ${src_files})