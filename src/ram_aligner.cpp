#include <string>
#include <vector>

#include "ram_aligner.hpp"

#define MATCH    0
#define INSERT   1
#define DELETE   2
#define MISMATCH 3
#define HEAD    -1

enum AlignmentType { global, semi_global, local };

typedef struct {
    int cost;
    int parent;
} cell;

void init_row(cell **matrix, int columns, int gap) {
    for (int j = 1; j < columns; j++) {
        matrix[0][j].cost = j * gap;
        matrix[0][j].parent = HEAD;
    }
}

void init_column(cell **matrix, int rows, int gap) {
    for (int i = 0; i < rows; i++) {
        matrix[i][0].cost = i * gap;
        matrix[i][0].parent = HEAD;
    }
}

void init_matrix(cell **matrix, AlignmentType type,
                 int rows, int columns, int gap) {
    matrix[0][0].cost = 0;
    matrix[0][0].parent = HEAD;

    switch (type) {
        case global:
            init_row(matrix, columns, gap);
            init_column(matrix, rows, gap);
            break;

            //prefix-suffix
        case semi_global:
            init_row(matrix, columns, 0);
            init_column(matrix, rows, gap);
            break;

        case local:
            init_row(matrix, columns, 0);
            init_column(matrix, rows, 0);
            break;
    }
}

int compare_string(cell **matrix, AlignmentType type,
                   const char *s, const char *t,
                   unsigned rows, unsigned columns,
                   int match, int mismatch, int gap,
                   unsigned *row, unsigned *col) {
    unsigned i, j, k;
    int options[3];
    int max_cost = 0;
    *row = 0;
    *col = 0;

    if (type == semi_global) {
        *row = 0;
        *col = columns - 1;
    }

    for (i = 1; i < rows; i++) {
        for (j = 1; j < columns; j++) {
            int m = (s[i - 1] == t[j - 1]) ? match : mismatch;

            options[MATCH] = matrix[i - 1][j - 1].cost + m;
            options[INSERT] = matrix[i][j - 1].cost + gap;
            options[DELETE] = matrix[i - 1][j].cost + gap;

            matrix[i][j].cost = options[MATCH];
            matrix[i][j].parent = (m == match) ? MATCH : MISMATCH;

            for (k = 1; k <= 2; k++) {
                if (options[k] > matrix[i][j].cost) {
                    matrix[i][j].cost = options[k];
                    matrix[i][j].parent = k;
                }
            }

            switch (type) {
                case local:
                    if (matrix[i][j].cost < 0) {
                        matrix[i][j].cost = 0;
                    }
                    if (max_cost < matrix[i][j].cost) {
                        max_cost = matrix[i][j].cost;
                        *row = i;
                        *col = j;
                    }
                    break;

                case semi_global:
                    if (j == columns - 1 && max_cost <= matrix[i][j].cost) {
                        max_cost = matrix[i][j].cost;
                        *row = i;
                        *col = j;
                    }
                    break;

                default:
                    break;
            }
        }
    }

    if (type == global) {
        *row = rows - 1;
        *col = columns - 1;
        return matrix[rows - 1][columns - 1].cost;

    } else {
        return max_cost;
    }
}

char get_letter(cell **matrix, unsigned row, unsigned col) {
    switch (matrix[row][col].parent) {
        case MATCH:
            return '=';
        case MISMATCH:
            return 'X';
        case INSERT:
            return 'I';
        case DELETE:
            return 'D';
        default:
            return ' ';
    }
}

void switch_cell(char letter, unsigned *row, unsigned *col) {
    switch (letter) {
        case 'I':
            (*col)--;
            break;
        case 'D':
            (*row)--;
            break;
        default:
            (*row)--;
            (*col)--;
            break;
    }
}

void make_cigar(cell **matrix, AlignmentType type,
                std::string &cigar, uint64_t &target_begin,
                unsigned row, unsigned col) {
    switch (type) {
        case local:
            while (matrix[row][col].cost != 0) {
                char letter = get_letter(matrix, row, col);
                cigar.push_back(letter);
                switch_cell(letter, &row, &col);
            }
            break;

        default:
            while (matrix[row][col].parent != HEAD) {
                char letter = get_letter(matrix, row, col);
                cigar.push_back(letter);
                switch_cell(letter, &row, &col);
            }

            while (row != 0) {
                row--;
                cigar.push_back('D');
            }

            while (col != 0) {
                col--;
                cigar.push_back('I');
            }
    }


}

int pairwise_alignment_connect(const char *query, uint64_t query_length,
                               const char *target, uint64_t target_length,
                               AlignmentType type,
                               int match, int mismatch, int gap,
                               std::string &cigar, uint64_t &target_begin, bool isCigar) {
    unsigned rows = query_length + 1;
    unsigned columns = target_length + 1;

    cell **matrix = new cell *[rows];
    for (unsigned i = 0; i < rows; i++) {
        matrix[i] = new cell[columns];
    }

    init_matrix(matrix, type, rows, columns, gap);

    unsigned row;
    unsigned col;

    int cost = compare_string(matrix, type, query, target, rows, columns, match, mismatch, gap, &row, &col);

    if (isCigar) {
        make_cigar(matrix, type, cigar, target_begin, row, col);
    }

    for (unsigned i = 0; i < rows; i++) {
        delete[] matrix[i];
    }

    delete[] matrix;

    return cost;
}

int pairwise_alignment(const char *query, uint64_t query_length,
                       const char *target, uint64_t target_length,
                       int match, int mismatch, int gap,
                       std::string &cigar, uint64_t &target_begin) {

    return pairwise_alignment_connect(query, query_length, target, target_length, local, match, mismatch, gap, cigar,
                                      target_begin, true);
}
