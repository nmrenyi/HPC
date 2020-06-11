// studentid: 2018011423
#include "util.h"
# define BLOCK_SIZE 32

// reference to https://docs.nvidia.com/cuda/cuda-c-programming-guide/#shared-memory
template<class T>
__global__ void my_gemm(T *A, T *B, T *C, int m, int n, int k, T alpha, T beta) {
    int blockRow = blockIdx.y;
    int blockCol = blockIdx.x;
	// find the start position for Csub
	T* Csub = &C[blockRow * BLOCK_SIZE * n + BLOCK_SIZE * blockCol];
	T Cvalue = 0.0;
    int row = threadIdx.y;
    int col = threadIdx.x;
	bool marginBlockDown = (blockRow == (m / BLOCK_SIZE));  // when m % BLOCK_SIZE != 0, identify the last bottom block
	bool marginBlockRight= (blockCol == (n / BLOCK_SIZE));  // when n % BLOCK_SIZE != 0, identify the last right block
	bool bottomMarginOK = (!marginBlockDown || row < (m % BLOCK_SIZE));  // when m % BLOCK_SIZE != 0, identify if the thread is within the valid block
	bool rightMarginOK = (!marginBlockRight || col < (n % BLOCK_SIZE));  // when n % BLOCK_SIZE != 0, identify if the thread is within the valid block
	bool okToCalculate = (bottomMarginOK && rightMarginOK);  // identify whether the current thread is valid for calculation(in case out of array range)

	// iterate over all needed subsections of A and B for calculating Csub
	// in this way, shared memory size can be fixed.
	__shared__ T As[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ T Bs[BLOCK_SIZE][BLOCK_SIZE];
	for (int i = 0; i < k / BLOCK_SIZE; i++) {
		T* Asub = &A[blockRow * BLOCK_SIZE * k + BLOCK_SIZE * i];
		T* Bsub = &B[i * BLOCK_SIZE * n + BLOCK_SIZE * blockCol];
		if (bottomMarginOK) {
			As[row][col] = Asub[row * k + col];
		}
		if (rightMarginOK) {
			Bs[row][col] = Bsub[row * n + col];
		}
		__syncthreads();
		if (okToCalculate) {
			for (int j = 0; j < BLOCK_SIZE; j++) {
				Cvalue += As[row][j] * Bs[j][col];
			}
		}
		__syncthreads();
	}

	// corner case: k is not divisible by BLOCK_SIZE
	int remaining = k % BLOCK_SIZE;
	if (remaining != 0) {
		T* Asub = &A[blockRow * BLOCK_SIZE * k + BLOCK_SIZE * (k / BLOCK_SIZE)];
		T* Bsub = &B[(k / BLOCK_SIZE) * BLOCK_SIZE * n + BLOCK_SIZE * blockCol];
		if (bottomMarginOK) {
			if (col < remaining)
				As[row][col] = Asub[row * k + col];
		}
		if (rightMarginOK) {
			if (row < remaining)
				Bs[row][col] = Bsub[row * n + col];
		}
		__syncthreads();
		if (okToCalculate) {
			for (int j = 0; j < remaining; j++) {
				Cvalue += As[row][j] * Bs[j][col];
			}
		}
	}

	// move value from register to global memory
	if (okToCalculate) {
		Csub[row * n + col] = Csub[row * n + col] * beta + Cvalue * alpha;
	}
}

template<class T>
double myGEMM(T* A, T* B, T* C, T alpha, T beta) {
	dim3 block(BLOCK_SIZE, BLOCK_SIZE);
	dim3 grid( (N + block.x - 1) / block.x, (M + block.y - 1) / block.y );
	my_gemm <<<grid, block>>>(A, B, C, M, N, K, alpha, beta);
	checkCudaErrors(cudaDeviceSynchronize());
	return 0.f;	
}
