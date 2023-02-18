#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_NODES 100
#define MAX_EDGES 1000
#define DAMPING_FACTOR 0.85
#define ERROR_THRESHOLD 0.0001

int n, m;
double PR[MAX_NODES];
int out_degree[MAX_NODES];
int graph[MAX_NODES][MAX_EDGES];

double calculate_new_PR(int node) {
  double new_PR = (1 - DAMPING_FACTOR);
  for (int j = 0; j < n; j++) {
    if (graph[j][node]) {
      new_PR += DAMPING_FACTOR * PR[j] / out_degree[j];
    }
  }
  return new_PR;
}

void PageRank() {
  int iteration = 0;
  double diff = 1;
  while (diff > ERROR_THRESHOLD) {
    diff = 0;
    for (int i = 0; i < n; i++) {
      double new_PR = calculate_new_PR(i);
      diff += fabs(PR[i] - new_PR);
      PR[i] = new_PR;
    }
    iteration++;
  }
  printf("PageRank converged after %d iterations.\n", iteration);
}

int main() {
  FILE *fp;
  fp = fopen("input.txt", "r");
  if (fp == NULL) {
    printf("Error opening file!\n");
    exit(1);
  }

  // Read the graph from file
  fscanf(fp, "%d%d", &n, &m);
  for (int i = 0; i < m; i++) {
    int from, to;
    fscanf(fp, "%d%d", &from, &to);
    graph[from][to] = 1;
    out_degree[from]++;
  }
  fclose(fp);

  // Initialize PageRank values
  for (int i = 0; i < n; i++) {
    PR[i] = 1.0 / n;
  }

  // Calculate PageRank
  PageRank();

  // Print the final PageRank values
  printf("Final PageRank values:\n");
  for (int i = 0; i < n; i++) {
    printf("PR[%d] = %lf\n", i, PR[i]);
  }
  return 0;
}