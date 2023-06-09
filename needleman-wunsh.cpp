#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <chrono>
#include <thread>
#include <atomic>
#include <cmath>
#include <mutex>
#include <map>

#define MAX_THREADS 5

void Needleman_Wusch_Thread(int** score_matrix, const std::string& S1, const std::string& S2,
                            int miss, int match, int gap,
                            char s1, char s2, int i, int j, int t, int l, int tl) {
    std::mutex lock;
    int score = (s1 != s2) ? miss : match;
    int max_score = std::max(std::max(t + gap, l + gap), tl + score);
    lock.lock();
    score_matrix[i][j] = max_score;
    lock.unlock();
}

void Needleman_Wusch(int** score_matrix, const std::vector<std::pair<int, int>>& matrix_points,
                     const std::string& S1, const std::string& S2,
                     int miss, int match, int gap) {
    for (const auto& mat_index : matrix_points) {
        char s1 = S1[mat_index.second - 1];
        char s2 = S2[mat_index.first - 1];
        Needleman_Wusch_Thread(score_matrix, S1, S2, miss, match, gap, s1, s2, mat_index.first, mat_index.second,
                               score_matrix[mat_index.first - 1][mat_index.second],
                               score_matrix[mat_index.first][mat_index.second - 1],
                               score_matrix[mat_index.first - 1][mat_index.second - 1]);
    }
}

std::vector<std::pair<char, char>> get_alignment(int** score_matrix, const std::string& S1, const std::string& S2,
                                                 int miss, int match, int gap) {
    std::vector<std::pair<char, char>> alignment;
    int i = S2.length();
    int j = S1.length();

    while (i > 0 && j > 0) {
        int score = score_matrix[i][j];
        int scoreDiag = score_matrix[i - 1][j - 1];
        int scoreLeft = score_matrix[i][j - 1];
        int scoreUp = score_matrix[i - 1][j];

        if (score == scoreDiag + ((S2[i - 1] != S1[j - 1]) ? miss : match)) {
            alignment.emplace_back(S2[i - 1], S1[j - 1]);
            i--;
            j--;
        }
        else if (score == scoreLeft + gap) {
            alignment.emplace_back('-', S1[j - 1]);
            j--;
        }
        else if (score == scoreUp + gap) {
            alignment.emplace_back(S2[i - 1], '-');
            i--;
        }
    }

    while (i > 0) {
        alignment.emplace_back(S2[i - 1], '-');
        i--;
    }
    while (j > 0) {
        alignment.emplace_back('-', S1[j - 1]);
        j--;
    }

    std::reverse(alignment.begin(), alignment.end());
    return alignment;
}

int main() {

    std::string S1 = "CCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCC";
    std::string S2 = "CCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCCTGTCC";
    auto start = std::chrono::high_resolution_clock::now();

    int miss = -1;
    int match = 1;
    int gap = -1;

    int m = S1.length() + 1;
    int n = S2.length() + 1;
    int** score_matrix = new int*[n];
    for (int i = 0; i < n; i++) {
        score_matrix[i] = new int[m];
        score_matrix[i][0] = i * gap;
    }
    for (int i = 0; i < m; i++) {
        score_matrix[0][i] = i * gap;
    }

    std::vector<std::vector<std::pair<int, int>>> matrix_points;
    for (int i = 2; i < m + 1; i++) {
        std::vector<std::pair<int, int>> curr_indices;
        for (int j = 1; j < i; j++) {
            if (j < n) {
                curr_indices.emplace_back(j, i - j);
            }
        }
        matrix_points.push_back(curr_indices);
    }
    for (int i = m + 1; i < m + n - 1; i++) {
        std::vector<std::pair<int, int>> curr_indices;
        for (int j = m - 1; j > i - n; j--) {
            curr_indices.emplace_back(i - j, j);
        }
        matrix_points.push_back(curr_indices);
    }

    for (int d = 0; d < matrix_points.size(); d++) {
        int diagonal = matrix_points[d].size();
        int num_threads = std::min(MAX_THREADS, diagonal);

        std::vector<std::thread> threads(num_threads);
        int block_size = diagonal / num_threads;

        std::vector<std::pair<int, int>> points1;
        for (int th = 0; th < num_threads - 1; th++) {
            for (int j = block_size * th; j < block_size * (th + 1); j++) {
                points1.push_back(matrix_points[d][j]);
            }
            threads[th] = std::thread(&Needleman_Wusch, score_matrix, std::ref(matrix_points[d]), std::cref(S1), std::cref(S2),
                                      miss, match, gap);
        }
        std::vector<std::pair<int, int>> points2;
        for (int j = diagonal - num_threads * block_size; j < diagonal; j++) {
            points2.push_back(matrix_points[d][j]);
        }

        threads[num_threads - 1] = std::thread(&Needleman_Wusch, score_matrix, std::ref(points2), std::cref(S1), std::cref(S2),
                                               miss, match, gap);

        for (int i = 0; i < num_threads; i++) {
            threads[i].join();
        }
    }

    auto elapsed = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = elapsed - start;
    std::cout << "Execution time: " << duration.count() << " seconds." << std::endl;

    std::vector<std::pair<char, char>> alignment = get_alignment(score_matrix, S1, S2, miss, match, gap);

    std::cout << "Optimal alignment:" << std::endl;
    for (const auto& pair : alignment) {
       // std::cout << pair.first << " " << pair.second << std::endl;
    }

    for (int i = 0; i < n; i++) {
        delete[] score_matrix[i];
    }
    delete[] score_matrix;

    return 0;
}