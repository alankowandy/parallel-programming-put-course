#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int isPrimeNumber(int number)
{
    int sqrtNumber = sqrt(number) + 1;
    for (int i = 2; i < sqrtNumber; ++i)
    {
        if (number % i == 0)
        {
            return 1;
        }
    }
    return 0;
}

void sequentiallyFindPrimeNumbers(int start, int MAX, int *array)
{
    for (int i = start; i <= MAX; i++)
    {
        array[i] = 1;
    }

    for (int i = start; i <= MAX; ++i)
    {
        if (isPrimeNumber(i) == 1)
        {
            array[i] = 0;
        }
    }
}

void parallelFindPrimeNumbers(int start, int MAX, int *array)
{
    #pragma omp parallel num_threads(12)
    {
        #pragma omp for
        for (int i = start; i <= MAX; i++)
        {
            array[i] = 1;
        }

        #pragma omp for schedule(dynamic)
        for (int i = start; i <= MAX; ++i)
        {
            if (isPrimeNumber(i) == 1)
            {
                array[i] = 0;
            }
        }
    }
}

void sequentialSieve(int* array, int max, int root) {
    for (int i = 2; i <= max; i++) {
        array[i] = 1;
    }

    for (int i = 2; i <= root; i++) {
        if (array[i] == 1) {
            for (int j = i * i; j <= max; j += i) {
                if (array[j] == 1)
                array[j] = 0;
            }
        }
    }
}

void domainSieve(int* array, int max, int root) {
    #pragma omp parallel num_threads(12)
    {
        #pragma omp for
        for (int i = 2; i <= max; i++) {
            array[i] = 1;
        }

        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();

        int start = 2 + thread_id;

        for (int i = start; i <= root; i += num_threads) {
            if (array[i] == 1) {
                for (int j = i * i; j <= max; j += i) {
                    array[j] = 0;
                }
            }
        }
    }
}

void functionSieve(int* array, int max, int root) {
    #pragma omp parallel num_threads(12)
    {
        #pragma omp for schedule(static)
        for (int i = 2; i <= max; i++) {
            array[i] = 1;
        }

        for (int i = 2; i <= root; i++) {
            if (array[i] == 1) {
                #pragma omp for schedule(dynamic)
                for (int j = i * i; j <= max; j += i) {
                    if (array[j] == 1)
                        array[j] = 0;
                }
            }
        }
    }
}

int compareArrays(int* arr1, int* arr2, int size) {
    for (int i = 2; i <= size; i++) {
        if (arr1[i] != arr2[i]) {
            return 0; // Arrays are different
        }
    }
    return 1; // Arrays are the same
}

void printArray(const int *array, int min, int max) {
    for (int i = min; i <= max; i++) {
        if (array[i] == 1) {
            printf("%d\t ", i);
        }
    }
    printf("\n");
}

int main() {
    int MAX = 100000000;
    int MIN = 2;
    int root = sqrt(MAX);

    double div_start, div_end, divp_start, divp_end;
    double seq_start, seq_end, domain_start, domain_end, func_start, func_end;

    int *div, *divp, *seq, *domain, *func;  // Array to mark composite numbers
    int div_count = 0, divp_count = 0, seq_count = 0, domain_count = 0, func_count = 0;  // Counter for composite numbers

    // Dynamic allocation of the arrays
    div = (int*)malloc((MAX + 1) * sizeof(int));
    divp = (int*)malloc((MAX + 1) * sizeof(int));
    seq = (int*)malloc((MAX + 1) * sizeof(int));
    func = (int*)malloc((MAX + 1) * sizeof(int));
    domain = (int*)malloc((MAX + 1) * sizeof(int));

    if (seq == NULL || func == NULL || domain == NULL || div == NULL || divp == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1; // Exit with an error code
    }

    // DZIELENIE SEKWENCYJNE

    div_start = omp_get_wtime();
    sequentiallyFindPrimeNumbers(MIN, MAX, div);
    div_end = omp_get_wtime();

    // DZIELENIE RÓWNOLEGŁE
    
    divp_start = omp_get_wtime();
    parallelFindPrimeNumbers(MIN, MAX, divp);
    divp_end = omp_get_wtime();
    
    //SEKWENCYJNIE

    seq_start = omp_get_wtime();
    sequentialSieve(seq, MAX, root);
    seq_end = omp_get_wtime();

    //RÓWNOLEGŁE FUNKCYJNE

    func_start = omp_get_wtime();
    functionSieve(func, MAX, root);
    func_end = omp_get_wtime();
    
    //RÓWNOLEGŁE DOMENOWO

    domain_start = omp_get_wtime();
    domainSieve(domain, MAX, root);
    domain_end = omp_get_wtime();

    // Wypisanie wyników i czasów

    for (int i = MIN; i <= MAX; i++) {
        div_count += (div[i] == 1);
        divp_count += (divp[i] == 1);
        seq_count += (seq[i] == 1);
        domain_count += (domain[i] == 1);
        func_count += (func[i] == 1);
    }
    
    printf("Instancja [%d, %d]", MIN, MAX);

    // Wypisanie wyników i czasów
    printf("\n| Method               | Wallclock Time (seconds) | Count of Prime Numbers |\n");
    printf("|----------------------|--------------------------|------------------------|\n");
    printf("| Sequential Division  | %-24.6f | %-22d |\n", div_end - div_start, div_count);
    printf("| Parallel Division    | %-24.6f | %-22d |\n", divp_end - divp_start, divp_count);
    printf("| Sequential Sieve     | %-24.6f | %-22d |\n", seq_end - seq_start, seq_count);
    printf("| Parallel Functional  | %-24.6f | %-22d |\n", func_end - func_start, func_count);
    printf("| Parallel Domain      | %-24.6f | %-22d |\n", domain_end - domain_start, domain_count);

    printf("\n");
    // Przyspieszenie
    printf("| Przyspieszenie (Acceleration) - Parallel Division:   | %-16.3f |\n", (div_end - div_start) / (divp_end - divp_start));
    printf("| Przyspieszenie (Acceleration) - Parallel Functional: | %-16.3f |\n", (seq_end - seq_start) / (func_end - func_start));
    printf("| Przyspieszenie (Acceleration) - Parallel Domain:     | %-16.3f |\n\n", (seq_end - seq_start) / (domain_end - domain_start));

    // Prędkość przetwarzania
    printf("| Speed - Sequential Division:       | %-28.2e/s |\n", (MAX - MIN) / (div_end - div_start));
    printf("| Speed - Parallel Division:         | %-28.2e/s |\n", (MAX - MIN) / (divp_end - divp_start));
    printf("| Speed - Sequential Sieve:          | %-28.2e/s |\n", (MAX - MIN) / (seq_end - seq_start));
    printf("| Speed - Parallel Functional:       | %-28.2e/s |\n", (MAX - MIN) / (func_end - func_start));
    printf("| Speed - Parallel Domain:           | %-28.2e/s |\n\n", (MAX - MIN) / (domain_end - domain_start));

    // Efektywność
    int num_threads = 12;
    double efficiency_divp = (div_end - div_start) / (num_threads * (divp_end - divp_start));
    double efficiency_func = (seq_end - seq_start) / (num_threads * (func_end - func_start));
    double efficiency_domain = (seq_end - seq_start) / (num_threads * (domain_end - domain_start));

    printf("| Efficiency - Parallel Division:                      | %-23.3f%% |\n", efficiency_divp * 100);
    printf("| Efficiency - Parallel Functional:                    | %-23.3f%% |\n", efficiency_func * 100);
    printf("| Efficiency - Parallel Domain:                        | %-23.3f%% |\n\n", efficiency_domain * 100);

    // Porównanie tabeli danych
    if (compareArrays(seq, func, MAX) && compareArrays(seq, domain, MAX) && compareArrays(seq, divp, MAX)) {
        printf("Tablice danych są takie same.\n");
    } else {
        printf("Tablice danych są różne.\n");
    }

    // printArray(div, MIN, MAX);
    // printArray(divp, MIN, MAX);
    // printArray(seq, MIN, MAX);
    // printArray(func, MIN, MAX);
    // printArray(domain, MIN, MAX);

    // Zwolnienie pamięci
    free(div);
    free(divp);
    free(seq);
    free(func);
    free(domain);

    return 0;
}