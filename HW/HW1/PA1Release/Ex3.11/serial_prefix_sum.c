#include <stdio.h>
#include <stdlib.h>

const int TOTAL_SIZE = 40;
int main() {
    int arr[TOTAL_SIZE];
    int i = 0;
    for (i = 0; i < TOTAL_SIZE; i++) {
        scanf("%d", &arr[i]);
    }
    int result = 0;
    for (i = 0; i < TOTAL_SIZE; i++) {
        result += arr[i];
        printf("%d ", result);
    }
    printf("\n");
    return 0;
}
