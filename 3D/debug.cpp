#include "debug.h"

// print three values of the array (the array should be of size n)
void print_fl_array(float* arr, int n, string label="") {
    label = (label.empty()) ? "" : "[" + label + "] ";

    int idx1 = (int) n / 4, idx2 = (int) n / 2, idx3 = (int) 3 * n / 4;
    float val1 = arr[idx1], val2 = arr[idx2], val3 = arr[idx3];
    cout << "---" << endl;
    cout << label << "idx " << idx1 << ": " << val1 << "; idx "
            << idx2 << ": " << val2 << "; idx " << idx3 << ": " << val3 << endl;
}

void print_fl_array_perc(float* arr, int n, float k, string label="") {
    label = (label.empty()) ? "" : "[" + label + "] ";

    for (int i = 0; i < n; i += (int) 1 / k) {
        cout << arr[i] << " ";
    }
    cout << endl << endl;
}
