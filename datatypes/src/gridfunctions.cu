
#include "gridfunctions.cuh"

#ifdef USING_CUDA

#include <cuda_runtime.h>

template <size_t D>
using scalar_arr_t = scalar_t[D];

template <size_t D>
using complex_arr_t = complex_t[D];

__global__ void fillDataKernel1d(const scalar_t* src, scalar_t* dest,
                                 int srcLength, int destLength) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;

  if (x < srcLength) {
    int offsetX = (destLength - srcLength) / 2;
    int destIndex = x + offsetX;
    dest[destIndex] = src[x];
  }
}
__global__ void fillDataKernel2d(const scalar_t* src, scalar_t* dest,
                                 int srcWidth, int srcHeight, int destWidth,
                                 int destHeight) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;

  if (x < srcWidth && y < srcHeight) {
    int offsetX = (destWidth - srcWidth) / 2;
    int offsetY = (destHeight - srcHeight) / 2;
    int destIndex = (y + offsetY) * destWidth + (x + offsetX);
    int srcIndex = y * srcWidth + x;
    dest[destIndex] = src[srcIndex];
  }
}

__global__ void fillDataKernel3d(const scalar_t* src, scalar_t* dest,
                                 int srcWidth, int srcHeight, int srcDepth,
                                 int destWidth, int destHeight, int destDepth) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z;

  if (x < srcWidth && y < srcHeight && z < srcDepth) {
    int offsetX = (destWidth - srcWidth) / 2;
    int offsetY = (destHeight - srcHeight) / 2;
    int offsetZ = (destDepth - srcDepth) / 2;
    int destIndex = ((z + offsetZ) * destHeight + (y + offsetY)) * destWidth +
                    (x + offsetX);
    int srcIndex = (z * srcHeight + y) * srcWidth + x;
    dest[destIndex] = src[srcIndex];
  }
}

__global__ void fillDataKernel1d(const complex_t* src, complex_t* dest,
                                 int srcLength, int destLength) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;

  if (x < srcLength) {
    int offsetX = (destLength - srcLength) / 2;
    int destIndex = x + offsetX;
    dest[destIndex] = src[x];
  }
}
__global__ void fillDataKernel2d(const complex_t* src, complex_t* dest,
                                 int srcWidth, int srcHeight, int destWidth,
                                 int destHeight) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;

  if (x < srcWidth && y < srcHeight) {
    int offsetX = (destWidth - srcWidth) / 2;
    int offsetY = (destHeight - srcHeight) / 2;
    int destIndex = (y + offsetY) * destWidth + (x + offsetX);
    int srcIndex = y * srcWidth + x;
    dest[destIndex] = src[srcIndex];
  }
}

__global__ void fillDataKernel3d(const complex_t* src, complex_t* dest,
                                 int srcWidth, int srcHeight, int srcDepth,
                                 int destWidth, int destHeight, int destDepth) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z;

  if (x < srcWidth && y < srcHeight && z < srcDepth) {
    int offsetX = (destWidth - srcWidth) / 2;
    int offsetY = (destHeight - srcHeight) / 2;
    int offsetZ = (destDepth - srcDepth) / 2;
    int destIndex = ((z + offsetZ) * destHeight + (y + offsetY)) * destWidth +
                    (x + offsetX);
    int srcIndex = (z * srcHeight + y) * srcWidth + x;
    dest[destIndex] = src[srcIndex];
  }
}

__global__ void fillDataKernel1d(const any_vector_t<scalar_t, 1>* src,
                                 scalar_t* (&dest)[1], int srcLength,
                                 int destLength) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;

  if (x < srcLength) {
    int offsetX = (destLength - srcLength) / 2;
    int destIndex = x + offsetX;
    dest[destIndex][0] = src[0][x];
  }
}
__global__ void fillDataKernel2d(const any_vector_t<scalar_t, 2>* src,
                                 scalar_t* (&dest)[2], int srcWidth,
                                 int srcHeight, int destWidth, int destHeight) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;

  if (x < srcWidth && y < srcHeight) {
    int offsetX = (destWidth - srcWidth) / 2;
    int offsetY = (destHeight - srcHeight) / 2;
    int destIndex = (y + offsetY) * destWidth + (x + offsetX);
    int srcIndex = y * srcWidth + x;
    dest[destIndex][0] = src[0][srcIndex];
    dest[destIndex][1] = src[1][srcIndex];
  }
}

__global__ void fillDataKernel3d(const any_vector_t<scalar_t, 3>* src,
                                 scalar_t* (&dest)[3], int srcWidth,
                                 int srcHeight, int srcDepth, int destWidth,
                                 int destHeight, int destDepth) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z;

  if (x < srcWidth && y < srcHeight && z < srcDepth) {
    int offsetX = (destWidth - srcWidth) / 2;
    int offsetY = (destHeight - srcHeight) / 2;
    int offsetZ = (destDepth - srcDepth) / 2;
    int destIndex = ((z + offsetZ) * destHeight + (y + offsetY)) * destWidth +
                    (x + offsetX);
    int srcIndex = (z * srcHeight + y) * srcWidth + x;
    dest[destIndex][0] = src[0][srcIndex];
    dest[destIndex][1] = src[1][srcIndex];
    dest[destIndex][2] = src[2][srcIndex];
  }
}

__global__ void fillDataKernel1d(const any_vector_t<complex_t, 1>* src,
                                 complex_t* (&dest)[1], int srcLength,
                                 int destLength) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;

  if (x < srcLength) {
    int offsetX = (destLength - srcLength) / 2;
    int destIndex = x + offsetX;
    dest[destIndex][0] = src[0][x];
  }
}
__global__ void fillDataKernel2d(const any_vector_t<complex_t, 2>* src,
                                 complex_t* (&dest)[2], int srcWidth,
                                 int srcHeight, int destWidth, int destHeight) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;

  if (x < srcWidth && y < srcHeight) {
    int offsetX = (destWidth - srcWidth) / 2;
    int offsetY = (destHeight - srcHeight) / 2;
    int destIndex = (y + offsetY) * destWidth + (x + offsetX);
    int srcIndex = y * srcWidth + x;
    dest[destIndex][0] = src[0][srcIndex];
    dest[destIndex][1] = src[1][srcIndex];
  }
}

__global__ void fillDataKernel3d(const any_vector_t<complex_t, 3>* src,
                                 complex_t* (&dest)[3], int srcWidth,
                                 int srcHeight, int srcDepth, int destWidth,
                                 int destHeight, int destDepth) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z;

  if (x < srcWidth && y < srcHeight && z < srcDepth) {
    int offsetX = (destWidth - srcWidth) / 2;
    int offsetY = (destHeight - srcHeight) / 2;
    int offsetZ = (destDepth - srcDepth) / 2;
    int destIndex = ((z + offsetZ) * destHeight + (y + offsetY)) * destWidth +
                    (x + offsetX);
    int srcIndex = (z * srcHeight + y) * srcWidth + x;
    dest[destIndex][0] = src[0][srcIndex];
    dest[destIndex][1] = src[1][srcIndex];
    dest[destIndex][2] = src[2][srcIndex];
  }
}

__global__ void centerDataKernel1d(const scalar_t* src, scalar_t* dest,
                                   int srcLength, int destLength) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;

  if (x < destLength) {
    int offsetX = (srcLength - destLength) / 2;
    int srcIndex = x + offsetX;
    dest[x] = src[srcIndex];
  }
}
__global__ void centerDataKernel2d(const scalar_t* src, scalar_t* dest,
                                   int srcWidth, int srcHeight, int destWidth,
                                   int destHeight) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;

  if (x < destWidth && y < destHeight) {
    int offsetX = (srcWidth - destWidth) / 2;
    int offsetY = (srcHeight - destHeight) / 2;
    int srcIndex = (y + offsetY) * srcWidth + (x + offsetX);
    int destIndex = y * destWidth + x;
    dest[destIndex] = src[srcIndex];
  }
}

__global__ void centerDataKernel3d(const scalar_t* src, scalar_t* dest,
                                   int srcWidth, int srcHeight, int srcDepth,
                                   int destWidth, int destHeight,
                                   int destDepth) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z;

  if (x < destWidth && y < destHeight && z < destDepth) {
    int offsetX = (srcWidth - destWidth) / 2;
    int offsetY = (srcHeight - destHeight) / 2;
    int offsetZ = (srcDepth - destDepth) / 2;
    int srcIndex =
        ((z + offsetZ) * srcHeight + (y + offsetY)) * srcWidth + (x + offsetX);
    int destIndex = (z * destHeight + y) * destWidth + x;
    dest[destIndex] = src[srcIndex];
  }
}

__global__ void centerDataKernel1d(const complex_t* src, complex_t* dest,
                                   int srcLength, int destLength) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;

  if (x < destLength) {
    int offsetX = (srcLength - destLength) / 2;
    int srcIndex = x + offsetX;
    dest[x] = src[srcIndex];
  }
}
__global__ void centerDataKernel2d(const complex_t* src, complex_t* dest,
                                   int srcWidth, int srcHeight, int destWidth,
                                   int destHeight) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;

  if (x < destWidth && y < destHeight) {
    int offsetX = (srcWidth - destWidth) / 2;
    int offsetY = (srcHeight - destHeight) / 2;
    int srcIndex = (y + offsetY) * srcWidth + (x + offsetX);
    int destIndex = y * destWidth + x;
    dest[destIndex] = src[srcIndex];
  }
}

__global__ void centerDataKernel3d(const complex_t* src, complex_t* dest,
                                   int srcWidth, int srcHeight, int srcDepth,
                                   int destWidth, int destHeight,
                                   int destDepth) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z;

  if (x < destWidth && y < destHeight && z < destDepth) {
    int offsetX = (srcWidth - destWidth) / 2;
    int offsetY = (srcHeight - destHeight) / 2;
    int offsetZ = (srcDepth - destDepth) / 2;
    int srcIndex =
        ((z + offsetZ) * srcHeight + (y + offsetY)) * srcWidth + (x + offsetX);
    int destIndex = (z * destHeight + y) * destWidth + x;
    dest[destIndex] = src[srcIndex];
  }
}

__global__ void centerDataKernel1d(scalar_t* src0,
                                   any_vector_t<scalar_t, 1>* dest,
                                   int srcLength, int destLength) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;

  if (x < destLength) {
    int offsetX = (srcLength - destLength) / 2;
    int srcIndex = x + offsetX;
    dest[x][0] = src0[srcIndex];
  }
}

__global__ void centerDataKernel2d(const scalar_t* src0, const scalar_t* src1,
                                   any_vector_t<scalar_t, 2>* dest,
                                   int srcWidth, int srcHeight, int destWidth,
                                   int destHeight) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;

  if (x < destWidth && y < destHeight) {
    int offsetX = (srcWidth - destWidth) / 2;
    int offsetY = (srcHeight - destHeight) / 2;
    int srcIndex = (y + offsetY) * srcWidth + (x + offsetX);
    int destIndex = y * destWidth + x;
    dest[destIndex][0] = src0[srcIndex];
    dest[destIndex][1] = src1[srcIndex];
  }
}

__global__ void centerDataKernel3d(const scalar_t* src0, const scalar_t* src1,
                                   const scalar_t* src2,
                                   any_vector_t<scalar_t, 3>* dest,
                                   int srcWidth, int srcHeight, int srcDepth,
                                   int destWidth, int destHeight,
                                   int destDepth) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z;

  if (x < destWidth && y < destHeight && z < destDepth) {
    int offsetX = (srcWidth - destWidth) / 2;
    int offsetY = (srcHeight - destHeight) / 2;
    int offsetZ = (srcDepth - destDepth) / 2;
    int srcIndex =
        ((z + offsetZ) * srcHeight + (y + offsetY)) * srcWidth + (x + offsetX);
    int destIndex = (z * destHeight + y) * destWidth + x;
    dest[destIndex][0] = src0[srcIndex];
    dest[destIndex][1] = src1[srcIndex];
    dest[destIndex][2] = src2[srcIndex];
  }
}

__global__ void centerDataKernel1d(const complex_t* src0,
                                   any_vector_t<complex_t, 1>* dest,
                                   int srcLength, int destLength) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;

  if (x < destLength) {
    int offsetX = (srcLength - destLength) / 2;
    int srcIndex = x + offsetX;
    dest[x][0] = src0[srcIndex];
  }
}
__global__ void centerDataKernel2d(const complex_t* src0, const complex_t* src1,
                                   any_vector_t<complex_t, 2>* dest,
                                   int srcWidth, int srcHeight, int destWidth,
                                   int destHeight) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;

  if (x < destWidth && y < destHeight) {
    int offsetX = (srcWidth - destWidth) / 2;
    int offsetY = (srcHeight - destHeight) / 2;
    int srcIndex = (y + offsetY) * srcWidth + (x + offsetX);
    int destIndex = y * destWidth + x;
    dest[destIndex][0] = src0[srcIndex];
    dest[destIndex][1] = src1[srcIndex];
  }
}

__global__ void centerDataKernel3d(const complex_t* src0, const complex_t* src1,
                                   const complex_t* src2,
                                   any_vector_t<complex_t, 3>* dest,
                                   int srcWidth, int srcHeight, int srcDepth,
                                   int destWidth, int destHeight,
                                   int destDepth) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z;

  if (x < destWidth && y < destHeight && z < destDepth) {
    int offsetX = (srcWidth - destWidth) / 2;
    int offsetY = (srcHeight - destHeight) / 2;
    int offsetZ = (srcDepth - destDepth) / 2;
    int srcIndex =
        ((z + offsetZ) * srcHeight + (y + offsetY)) * srcWidth + (x + offsetX);
    int destIndex = (z * destHeight + y) * destWidth + x;
    dest[destIndex][0] = src0[srcIndex];
    dest[destIndex][1] = src1[srcIndex];
    dest[destIndex][2] = src2[srcIndex];
  }
}

void grid::fill_interior_cuda_1d(const scalar_t* srcHost, scalar_t* destDevice,
                                 int srcLength, int destLength) {
  scalar_t* srcDevice;
  CHECK_CUDA_ERROR(cudaMalloc(&srcDevice, srcLength * sizeof(scalar_t)));
  CHECK_CUDA_ERROR(cudaMemcpy(srcDevice, srcHost, srcLength * sizeof(scalar_t),
                              cudaMemcpyHostToDevice));

  int threadsPerBlock = 256;
  int blocksPerGrid = (srcLength + threadsPerBlock - 1) / threadsPerBlock;
  fillDataKernel1d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcLength, destLength);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(srcDevice));
}

void grid::fill_interior_cuda_2d(const scalar_t* srcHost, scalar_t* destDevice,
                                 int srcWidth, int srcHeight, int destWidth,
                                 int destHeight) {
  scalar_t* srcDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&srcDevice, srcWidth * srcHeight * sizeof(scalar_t)));
  CHECK_CUDA_ERROR(cudaMemcpy(srcDevice, srcHost,
                              srcWidth * srcHeight * sizeof(scalar_t),
                              cudaMemcpyHostToDevice));

  dim3 threadsPerBlock(16, 16);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y);
  fillDataKernel2d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcWidth, srcHeight, destWidth, destHeight);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(srcDevice));
}

void grid::fill_interior_cuda_3d(const scalar_t* srcHost, scalar_t* destDevice,
                                 int srcWidth, int srcHeight, int srcDepth,
                                 int destWidth, int destHeight, int destDepth) {
  scalar_t* srcDevice;
  CHECK_CUDA_ERROR(cudaMalloc(
      &srcDevice, srcWidth * srcHeight * srcDepth * sizeof(scalar_t)));
  CHECK_CUDA_ERROR(cudaMemcpy(
      srcDevice, srcHost, srcWidth * srcHeight * srcDepth * sizeof(scalar_t),
      cudaMemcpyHostToDevice));

  dim3 threadsPerBlock(8, 8, 8);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (srcDepth + threadsPerBlock.z - 1) / threadsPerBlock.z);

  // Launch the kernel
  fillDataKernel3d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcWidth, srcHeight, srcDepth, destWidth,
      destHeight, destDepth);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(srcDevice));
}

void grid::copy_interior_cuda_1d(const scalar_t* srcDevice, scalar_t* destHost,
                                 int srcLength, int destLength) {
  scalar_t* destDevice;
  CHECK_CUDA_ERROR(cudaMalloc(&destDevice, destLength * sizeof(scalar_t)));

  int threadsPerBlock = 256;
  int blocksPerGrid = (srcLength + threadsPerBlock - 1) / threadsPerBlock;
  centerDataKernel1d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcLength, destLength);

  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaMemcpy(destHost, destDevice,
                              destLength * sizeof(scalar_t),
                              cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(destDevice));
}

void grid::copy_interior_cuda_2d(const scalar_t* srcDevice, scalar_t* destHost,
                                 int srcWidth, int srcHeight, int destWidth,
                                 int destHeight) {
  scalar_t* destDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&destDevice, destWidth * destHeight * sizeof(scalar_t)));

  dim3 threadsPerBlock(16, 16);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y);
  centerDataKernel2d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcWidth, srcHeight, destWidth, destHeight);

  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaMemcpy(destHost, destDevice,
                              destWidth * destHeight * sizeof(scalar_t),
                              cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(destDevice));
}

void grid::copy_interior_cuda_3d(const scalar_t* srcDevice, scalar_t* destHost,
                                 int srcWidth, int srcHeight, int srcDepth,
                                 int destWidth, int destHeight, int destDepth) {
  scalar_t* destDevice;
  CHECK_CUDA_ERROR(cudaMalloc(
      &destDevice, srcWidth * srcHeight * srcDepth * sizeof(scalar_t)));

  dim3 threadsPerBlock(8, 8, 8);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (srcDepth + threadsPerBlock.z - 1) / threadsPerBlock.z);

  // Launch the kernel
  centerDataKernel3d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcWidth, srcHeight, srcDepth, destWidth,
      destHeight, destDepth);

  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(
      cudaMemcpy(destHost, destDevice,
                 destWidth * destHeight * destDepth * sizeof(scalar_t),
                 cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(destDevice));
}

void grid::fill_interior_cuda_1d(const complex_t* srcHost,
                                 complex_t* destDevice, int srcLength,
                                 int destLength) {
  complex_t* srcDevice;
  CHECK_CUDA_ERROR(cudaMalloc(&srcDevice, srcLength * sizeof(complex_t)));
  CHECK_CUDA_ERROR(cudaMemcpy(srcDevice, srcHost, srcLength * sizeof(complex_t),
                              cudaMemcpyHostToDevice));

  int threadsPerBlock = 256;
  int blocksPerGrid = (srcLength + threadsPerBlock - 1) / threadsPerBlock;
  fillDataKernel1d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcLength, destLength);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(srcDevice));
}

void grid::fill_interior_cuda_2d(const complex_t* srcHost,
                                 complex_t* destDevice, int srcWidth,
                                 int srcHeight, int destWidth, int destHeight) {
  complex_t* srcDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&srcDevice, srcWidth * srcHeight * sizeof(complex_t)));
  CHECK_CUDA_ERROR(cudaMemcpy(srcDevice, srcHost,
                              srcWidth * srcHeight * sizeof(complex_t),
                              cudaMemcpyHostToDevice));

  dim3 threadsPerBlock(16, 16);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y);
  fillDataKernel2d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcWidth, srcHeight, destWidth, destHeight);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(srcDevice));
}

void grid::fill_interior_cuda_3d(const complex_t* srcHost,
                                 complex_t* destDevice, int srcWidth,
                                 int srcHeight, int srcDepth, int destWidth,
                                 int destHeight, int destDepth) {
  complex_t* srcDevice;
  CHECK_CUDA_ERROR(cudaMalloc(
      &srcDevice, srcWidth * srcHeight * srcDepth * sizeof(complex_t)));
  CHECK_CUDA_ERROR(cudaMemcpy(
      srcDevice, srcHost, srcWidth * srcHeight * srcDepth * sizeof(complex_t),
      cudaMemcpyHostToDevice));

  dim3 threadsPerBlock(8, 8, 8);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (srcDepth + threadsPerBlock.z - 1) / threadsPerBlock.z);

  // Launch the kernel
  fillDataKernel3d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcWidth, srcHeight, srcDepth, destWidth,
      destHeight, destDepth);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(srcDevice));
}

void grid::copy_interior_cuda_1d(const complex_t* srcDevice,
                                 complex_t* destHost, int srcLength,
                                 int destLength) {
  complex_t* destDevice;
  CHECK_CUDA_ERROR(cudaMalloc(&destDevice, destLength * sizeof(complex_t)));

  int threadsPerBlock = 256;
  int blocksPerGrid = (srcLength + threadsPerBlock - 1) / threadsPerBlock;
  centerDataKernel1d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcLength, destLength);

  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaMemcpy(destHost, destDevice,
                              destLength * sizeof(complex_t),
                              cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(destDevice));
}

void grid::copy_interior_cuda_2d(const complex_t* srcDevice,
                                 complex_t* destHost, int srcWidth,
                                 int srcHeight, int destWidth, int destHeight) {
  complex_t* destDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&destDevice, destWidth * destHeight * sizeof(complex_t)));

  dim3 threadsPerBlock(16, 16);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y);
  centerDataKernel2d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcWidth, srcHeight, destWidth, destHeight);

  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaMemcpy(destHost, destDevice,
                              destWidth * destHeight * sizeof(complex_t),
                              cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(destDevice));
}

void grid::copy_interior_cuda_3d(const complex_t* srcDevice,
                                 complex_t* destHost, int srcWidth,
                                 int srcHeight, int srcDepth, int destWidth,
                                 int destHeight, int destDepth) {
  complex_t* destDevice;
  CHECK_CUDA_ERROR(cudaMalloc(
      &destDevice, srcWidth * srcHeight * srcDepth * sizeof(complex_t)));

  dim3 threadsPerBlock(8, 8, 8);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (srcDepth + threadsPerBlock.z - 1) / threadsPerBlock.z);

  // Launch the kernel
  centerDataKernel3d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcWidth, srcHeight, srcDepth, destWidth,
      destHeight, destDepth);

  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(
      cudaMemcpy(destHost, destDevice,
                 destWidth * destHeight * destDepth * sizeof(complex_t),
                 cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(destDevice));
}

void grid::fill_interior_cuda_1d(const any_vector_t<scalar_t, 1>* srcHost,
                                 scalar_t* (&destDevice)[1], int srcLength,
                                 int destLength) {
  any_vector_t<scalar_t, 1>* srcDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&srcDevice, srcLength * sizeof(any_vector_t<scalar_t, 1>)));
  CHECK_CUDA_ERROR(cudaMemcpy(srcDevice, srcHost,
                              srcLength * sizeof(any_vector_t<scalar_t, 1>),
                              cudaMemcpyHostToDevice));

  int threadsPerBlock = 256;
  int blocksPerGrid = (srcLength + threadsPerBlock - 1) / threadsPerBlock;
  fillDataKernel1d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcLength, destLength);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(srcDevice));
}

void grid::fill_interior_cuda_2d(const any_vector_t<scalar_t, 2>* srcHost,
                                 scalar_t* (&destDevice)[2], int srcWidth,
                                 int srcHeight, int destWidth, int destHeight) {
  any_vector_t<scalar_t, 2>* srcDevice;
  CHECK_CUDA_ERROR(cudaMalloc(
      &srcDevice, srcWidth * srcHeight * sizeof(any_vector_t<scalar_t, 2>)));
  CHECK_CUDA_ERROR(
      cudaMemcpy(srcDevice, srcHost,
                 srcWidth * srcHeight * sizeof(any_vector_t<scalar_t, 2>),
                 cudaMemcpyHostToDevice));

  dim3 threadsPerBlock(16, 16);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y);
  fillDataKernel2d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcWidth, srcHeight, destWidth, destHeight);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(srcDevice));
}

void grid::fill_interior_cuda_3d(const any_vector_t<scalar_t, 3>* srcHost,
                                 scalar_t* (&destDevice)[3], int srcWidth,
                                 int srcHeight, int srcDepth, int destWidth,
                                 int destHeight, int destDepth) {
  any_vector_t<scalar_t, 3>* srcDevice;
  CHECK_CUDA_ERROR(cudaMalloc(
      &srcDevice, srcWidth * srcHeight * srcDepth * sizeof(scalar_t)));
  CHECK_CUDA_ERROR(cudaMemcpy(
      srcDevice, srcHost, srcWidth * srcHeight * srcDepth * sizeof(scalar_t),
      cudaMemcpyHostToDevice));

  dim3 threadsPerBlock(8, 8, 8);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (srcDepth + threadsPerBlock.z - 1) / threadsPerBlock.z);

  // Launch the kernel
  fillDataKernel3d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcWidth, srcHeight, srcDepth, destWidth,
      destHeight, destDepth);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(srcDevice));
}

void grid::copy_interior_cuda_1d(scalar_t* const (&srcDevice)[1],
                                 any_vector_t<scalar_t, 1>* destHost,
                                 int srcLength, int destLength) {
  any_vector_t<scalar_t, 1>* destDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&destDevice, destLength * sizeof(any_vector_t<scalar_t, 1>)));

  int threadsPerBlock = 256;
  int blocksPerGrid = (srcLength + threadsPerBlock - 1) / threadsPerBlock;
  centerDataKernel1d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice[0], destDevice, srcLength, destLength);

  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaMemcpy(destHost, destDevice,
                              destLength * sizeof(any_vector_t<scalar_t, 1>),
                              cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(destDevice));
}

void grid::copy_interior_cuda_2d(scalar_t* const (&srcDevice)[2],
                                 any_vector_t<scalar_t, 2>* destHost,
                                 int srcWidth, int srcHeight, int destWidth,
                                 int destHeight) {
  any_vector_t<scalar_t, 2>* destDevice;
  CHECK_CUDA_ERROR(cudaMalloc(
      &destDevice, destWidth * destHeight * sizeof(any_vector_t<scalar_t, 2>)));

  dim3 threadsPerBlock(16, 16);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y);
  centerDataKernel2d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice[0], srcDevice[1], destDevice, srcWidth, srcHeight, destWidth,
      destHeight);

  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaMemcpy(destHost, destDevice,
                              destWidth * destHeight * sizeof(scalar_arr_t<2>),
                              cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(destDevice));
}

void grid::copy_interior_cuda_3d(scalar_t* const (&srcDevice)[3],
                                 any_vector_t<scalar_t, 3>* destHost,
                                 int srcWidth, int srcHeight, int srcDepth,
                                 int destWidth, int destHeight, int destDepth) {
  any_vector_t<scalar_t, 3>* destDevice;
  CHECK_CUDA_ERROR(cudaMalloc(
      &destDevice, destWidth * destHeight * sizeof(any_vector_t<scalar_t, 3>)));

  dim3 threadsPerBlock(16, 16);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y);
  centerDataKernel3d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice[0], srcDevice[1], srcDevice[2], destDevice, srcWidth, srcHeight,
      srcDepth, destWidth, destHeight, destDepth);

  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(
      cudaMemcpy(destHost, destDevice,
                 destWidth * destHeight * sizeof(any_vector_t<scalar_t, 3>),
                 cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(destDevice));
}

void grid::fill_interior_cuda_1d(const any_vector_t<complex_t, 1>* srcHost,
                                 complex_t* (&destDevice)[1], int srcLength,
                                 int destLength) {
  any_vector_t<complex_t, 1>* srcDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&srcDevice, srcLength * sizeof(any_vector_t<complex_t, 1>)));
  CHECK_CUDA_ERROR(cudaMemcpy(srcDevice, srcHost,
                              srcLength * sizeof(any_vector_t<complex_t, 1>),
                              cudaMemcpyHostToDevice));

  int threadsPerBlock = 256;
  int blocksPerGrid = (srcLength + threadsPerBlock - 1) / threadsPerBlock;
  fillDataKernel1d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcLength, destLength);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(srcDevice));
}

void grid::fill_interior_cuda_2d(const any_vector_t<complex_t, 2>* srcHost,
                                 complex_t* (&destDevice)[2], int srcWidth,
                                 int srcHeight, int destWidth, int destHeight) {
  any_vector_t<complex_t, 2>* srcDevice;
  CHECK_CUDA_ERROR(cudaMalloc(
      &srcDevice, srcWidth * srcHeight * sizeof(any_vector_t<complex_t, 2>)));
  CHECK_CUDA_ERROR(
      cudaMemcpy(srcDevice, srcHost,
                 srcWidth * srcHeight * sizeof(any_vector_t<complex_t, 2>),
                 cudaMemcpyHostToDevice));

  dim3 threadsPerBlock(16, 16);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y);
  fillDataKernel2d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcWidth, srcHeight, destWidth, destHeight);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(srcDevice));
}

void grid::fill_interior_cuda_3d(const any_vector_t<complex_t, 3>* srcHost,
                                 complex_t* (&destDevice)[3], int srcWidth,
                                 int srcHeight, int srcDepth, int destWidth,
                                 int destHeight, int destDepth) {
  any_vector_t<complex_t, 3>* srcDevice;
  CHECK_CUDA_ERROR(cudaMalloc(
      &srcDevice, srcWidth * srcHeight * srcDepth * sizeof(complex_t)));
  CHECK_CUDA_ERROR(cudaMemcpy(
      srcDevice, srcHost, srcWidth * srcHeight * srcDepth * sizeof(complex_t),
      cudaMemcpyHostToDevice));

  dim3 threadsPerBlock(8, 8, 8);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y,
                     (srcDepth + threadsPerBlock.z - 1) / threadsPerBlock.z);

  // Launch the kernel
  fillDataKernel3d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice, destDevice, srcWidth, srcHeight, srcDepth, destWidth,
      destHeight, destDepth);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(srcDevice));
}

void grid::copy_interior_cuda_1d(complex_t* const (&srcDevice)[1],
                                 any_vector_t<complex_t, 1>* destHost,
                                 int srcLength, int destLength) {
  any_vector_t<complex_t, 1>* destDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&destDevice, destLength * sizeof(any_vector_t<complex_t, 1>)));

  int threadsPerBlock = 256;
  int blocksPerGrid = (srcLength + threadsPerBlock - 1) / threadsPerBlock;
  centerDataKernel1d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice[0], destDevice, srcLength, destLength);

  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaMemcpy(destHost, destDevice,
                              destLength * sizeof(any_vector_t<complex_t, 1>),
                              cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(destDevice));
}

void grid::copy_interior_cuda_2d(complex_t* const (&srcDevice)[2],
                                 any_vector_t<complex_t, 2>* destHost,
                                 int srcWidth, int srcHeight, int destWidth,
                                 int destHeight) {
  any_vector_t<complex_t, 2>* destDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&destDevice,
                 destWidth * destHeight * sizeof(any_vector_t<complex_t, 2>)));

  dim3 threadsPerBlock(16, 16);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y);
  centerDataKernel2d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice[0], srcDevice[1], destDevice, srcWidth, srcHeight, destWidth,
      destHeight);

  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(
      cudaMemcpy(destHost, destDevice,
                 destWidth * destHeight * sizeof(any_vector_t<complex_t, 2>),
                 cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(destDevice));
}

void grid::copy_interior_cuda_3d(complex_t* const (&srcDevice)[3],
                                 any_vector_t<complex_t, 3>* destHost,
                                 int srcWidth, int srcHeight, int srcDepth,
                                 int destWidth, int destHeight, int destDepth) {
  any_vector_t<complex_t, 3>* destDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&destDevice,
                 destWidth * destHeight * sizeof(any_vector_t<complex_t, 3>)));

  dim3 threadsPerBlock(16, 16);
  dim3 blocksPerGrid((srcWidth + threadsPerBlock.x - 1) / threadsPerBlock.x,
                     (srcHeight + threadsPerBlock.y - 1) / threadsPerBlock.y);
  centerDataKernel3d CUDA_KERNEL(blocksPerGrid, threadsPerBlock)(
      srcDevice[0], srcDevice[1], srcDevice[2], destDevice, srcWidth, srcHeight,
      srcDepth, destWidth, destHeight, destDepth);

  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(
      cudaMemcpy(destHost, destDevice,
                 destWidth * destHeight * sizeof(any_vector_t<complex_t, 3>),
                 cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(destDevice));
}

#endif