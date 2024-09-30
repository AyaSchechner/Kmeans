#include <stdio.h>
# include <stdlib.h>
#include <math.h>


struct cord{
    double value;
    struct cord *next;
};

struct vector{
    struct vector *next;
    struct cord *cords;
};


double epsilon = 0.001;
double euclidean_dis(double *, double *, int);
int check_centroids(double *, int);
void calc_new_cent(double **, int *, double **, double *, int , int *, int , int);
void adjust_centroids(double **, double **, int, int *, int *, int, int, double *);
int vectors_num(struct vector *);
int vectors_len(struct vector *);
double ** insert_vectors(double **, struct vector *, int, int);
double ** create_vectors_mat(int, int, struct vector *);
double ** create_centroids_mat(int, int, double **);
void first_assign_vectors(double **,  double **, int *, int *, int, int, int);
void assign_vectors(double **,  double **, int *, int *, int, int, int);
void print_centroids(double **, int, int);
void free_mat(double **, int);
void free_mem(double **, double **, int *, int *, double *, int, int);
void free_exit(char *, char *);
void free_LinkedList(struct vector *, int);
int k_means(int, int, struct vector *);



/* computes number of vectors in the input */
int vectors_num(struct vector * vectors){
    int N = 0;
    struct vector * curr = vectors;
    while (curr->next != NULL){
        N++;
        curr = curr->next;
    }
    return N;
}


/* computes the length of the vectors in the input */
int vectors_len(struct vector * vectors){
    int d = 0;
    struct vector * vec = vectors;
    struct cord * curr = vec->cords;
    while (curr != NULL){
        d++;
        curr = curr->next;
    }
    return d;
}


/* inserting vectors into the matrix from the linked list */
double ** insert_vectors(double ** matrix, struct vector * vectors_ll, int vec_num, int vec_len){
    int i;
    int j;
    struct vector *curr = vectors_ll;
    struct cord *curr1;
    for (i = 0; i < vec_num; ++i) {
        curr1 = curr->cords;
        for (j = 0; j < vec_len ; ++j) {
            matrix[i][j] = curr1->value;
            curr1 = curr1->next;
        }
        curr = curr->next;
    }
    return matrix;
}


/* creating matrix to store the vectors */
double ** create_vectors_mat(int vec_num, int vec_len, struct vector * vectors_ll){
    double **vectors_mat = (double **) malloc(vec_num * sizeof(double*));
    int i = 0;
    if(vectors_mat == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }
    while (i < vec_num){
        vectors_mat[i] = (double *) malloc(vec_len * sizeof(double));
        if(vectors_mat[i] == NULL){
            printf("An Error Has Occurred");
            free_mat(vectors_mat, i);
            return NULL;
        }
        i++;
    }
    vectors_mat = insert_vectors(vectors_mat, vectors_ll,vec_num, vec_len);
    return vectors_mat;
}


/* creating matrix to store centroids */
double ** create_centroids_mat(int k, int vec_len, double ** vectors_mat){
    int i;
    int j;
    double **centroids_mat = (double **) malloc(k * sizeof(double*));
    if(vectors_mat == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }
    for (i = 0; i < k; ++i) {
        centroids_mat[i] = (double *) malloc(vec_len * sizeof(double));
        if(vectors_mat == NULL){
            printf("An Error Has Occurred");
            free_mat(centroids_mat, i);
            return NULL;
        }
    }
    for (i = 0; i < k; ++i) {
        for (j = 0; j < vec_len; ++j) {
            centroids_mat[i][j] = vectors_mat[i][j];
        }
    }
    return centroids_mat;
}


/* computes Euclidean distance between two vectors */
double euclidean_dis(double * vec1, double * vec2, int vec_len) {
    int i;
    double dis = 0.0;
    for (i = 0; i < vec_len; i++) {
        dis += (((vec1[i]) - (vec2[i])) * ((vec1[i]) - (vec2[i])));
    }
    dis = sqrt(dis);
    return dis;
}


/* assign vectors to clusters for the first time */
void first_assign_vectors(double ** vectors_mat,  double ** centroids_mat, int * index, int * cluster_size, int vec_num, int vec_len, int k) {
    int i;
    int j;
    double dis_first_cent;
    double dis_j_cent;
    int curr_index;
    for (i = 0; i < vec_num; i++) {
        dis_first_cent = euclidean_dis(vectors_mat[i], centroids_mat[0], vec_len);
        curr_index = 0;
        for (j = 1; j < k ; j++) {
            dis_j_cent = euclidean_dis(vectors_mat[i], centroids_mat[j], vec_len);
            if (dis_first_cent > dis_j_cent) {
                dis_first_cent = dis_j_cent;
                curr_index = j;
            }
        }
        index[i] = curr_index;
        cluster_size[curr_index]++;
    }
}


/* computing each new centroid separately */
void calc_new_cent(double **vec_mat, int *index, double ** cent_mat, double *deltas, int cluster, int * cluster_size, int vec_num, int vec_len){
    int i;
    double delta;
    delta = 0.0;
    for (i = 0; i < vec_len; i++){
        int j;
        double curr_index;
        curr_index = 0.0;
        for (j = 0; j < vec_num; j++) {
            if(index[j] == cluster){
                curr_index += vec_mat[j][i];
            }
        }
        delta += ((cent_mat[cluster][i])-((curr_index)/(cluster_size[cluster])))*((cent_mat[cluster][i])-((curr_index)/(cluster_size[cluster])));

        cent_mat[cluster][i] = ((curr_index)/(cluster_size[cluster]));
    }
    deltas[cluster] = sqrt(delta);
}

/* computing all new centroids after adjusting vectors + computing the delta of each centroid after its change */
void adjust_centroids(double **vec_mat, double **cent_mat, int k, int *index, int * cluster_size, int vec_num, int vec_len, double *delta) {
    int i;
    for (i = 0; i < k; i++) {
        calc_new_cent(vec_mat, index, cent_mat, delta, i, cluster_size, vec_num, vec_len);
    }

}


/* checking whether all centroids deltas are lower then epsilon */
int check_centroids(double *centroids_delta, int k){
    int i;
    for (i = 0; i < k; i++) {
        if(centroids_delta[i] >= epsilon){
            return 0;
        }
    }
    return 1;
}


/* assigning vectors to the clusters (not the first time assigning!) */
void assign_vectors(double ** vectors_mat,  double ** centroids_mat, int * index, int * cluster_size, int vec_num, int vec_len, int k) {
    int i;
    int j;
    double dis_prev_cent;
    double dis_j_cent;
    for (i = 0; i < vec_num; i++) {
        dis_prev_cent = euclidean_dis(vectors_mat[i], centroids_mat[index[i]], vec_len);
        for (j = 0; j < k; j++) {
            if(j == index[i]){
                continue;
            }
            dis_j_cent = euclidean_dis(vectors_mat[i], centroids_mat[j], vec_len);
            if(dis_j_cent < dis_prev_cent){
                dis_prev_cent = dis_j_cent;
                cluster_size[index[i]]--;
                index[i] = j;
                cluster_size[j]++;
            }
        }
    }
}


/* print all centroids in case we finished */
void print_centroids(double ** centroids_mat, int vec_len, int k){
    int i;
    int j;
    int d;
    for (i = 0; i < k; i++) {
        d = vec_len;
        for (j = 0; j < vec_len; j++) {
            printf("%.4f", centroids_mat[i][j]);
            if (d > 1){
                printf(",");
                d--;
            }
            else{
                printf("\n");
            }
        }
    }
}

/* freeing memory */
void free_mat(double ** matrix, int rows_num) {
    int i;
    for (i = 0; i < rows_num; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

/* freeing memory */
void free_mem(double ** vectors_mat, double ** centroids_mat, int * index, int * cluster_size, double * centroids_deltas, int vec_num, int k){
    free_mat(centroids_mat, k);
    free_mat(vectors_mat, vec_num);
    free(cluster_size);
    free(centroids_deltas);
    free(index);
}

/* freeing memory */
void free_exit(char * cluster_num, char * iter){
    free(cluster_num);
    free(iter);
}


void free_LinkedList(struct vector *head_vec, int d) {
    struct vector *curr_vec = head_vec;
    struct vector *next_vec;
    int a = 0;

    while (a < d) {
        struct cord *curr_cord = curr_vec->cords;
        struct cord *next_cord;

        while (curr_cord != NULL) {
            next_cord = curr_cord->next;
            free(curr_cord);
            curr_cord = next_cord;
            a++;
        }

        next_vec = curr_vec->next;
        free(curr_vec);
        curr_vec = next_vec;
    }
    free(curr_vec);
}


/* k_means function */
int k_means(int k, int iter, struct vector * input_data){
    int vec_num = vectors_num(input_data);
    int vec_len = vectors_len(input_data);
    int c = vec_len*vec_num;
    int * index;
    int * cluster_size;
    double ** vectors_mat;
    double ** centroids_mat;
    int check;
    int iter_num;
    double * centroids_deltas;
    if(k >= vec_num){
        printf("Invalid number of clusters!");
        free_LinkedList(input_data, c);
        return 1;
    }
    vectors_mat = create_vectors_mat(vec_num, vec_len, input_data);
    if(vectors_mat == NULL){
        free_LinkedList(input_data, c);
        return 1;
    }
    free_LinkedList(input_data,c);
    centroids_mat = create_centroids_mat(k, vec_len, vectors_mat);
    if(centroids_mat == NULL){
        free_mat(vectors_mat,vec_num);
        return 1;
    }
    index = malloc(vec_num * sizeof(int));                            /* index[i] = j - vector i belongs to cluster j (size N) */
    if(index == NULL){
        printf("An Error Has Occurred");
        free_mat(vectors_mat,vec_num);
        free_mat(centroids_mat,k);
        return 1;
    }
    cluster_size = calloc(k , sizeof(int));    /* cluster_size[i] = j - there are j vectors in cluster i (size k) */
    centroids_deltas = malloc(k * sizeof(double ));               /* delta(prev cluster, current cluster) */
    if(centroids_deltas == NULL){
        printf("An Error Has Occurred");
        free_mat(vectors_mat,vec_num);
        free_mat(centroids_mat,k);
        free(index);
        return 1;
    }
    first_assign_vectors(vectors_mat, centroids_mat, index, cluster_size, vec_num, vec_len, k);
    adjust_centroids(vectors_mat, centroids_mat, k, index, cluster_size, vec_num, vec_len, centroids_deltas);
    check = check_centroids(centroids_deltas, k);
    if (check == 1){        /* if all deltas are smaller then epsilon */
        print_centroids(centroids_mat, vec_len, k);
        free_mem(vectors_mat, centroids_mat, index, cluster_size, centroids_deltas, vec_num, k);
        return 0;
    }
    iter_num = 1;
    while (iter_num < iter){
        assign_vectors(vectors_mat, centroids_mat, index, cluster_size, vec_num, vec_len, k);
        adjust_centroids(vectors_mat, centroids_mat, k, index, cluster_size, vec_num, vec_len, centroids_deltas);
        check = check_centroids(centroids_deltas, k);
        if (check == 1){
            print_centroids(centroids_mat, vec_len, k);
            free_mem(vectors_mat, centroids_mat, index, cluster_size, centroids_deltas, vec_num, k);
            return 0;
        }
        iter_num++;
    }
    print_centroids(centroids_mat, vec_len, k);
    free_mem(vectors_mat, centroids_mat, index, cluster_size, centroids_deltas, vec_num, k);
    return 0;
}


int main(int argc, char **argv){
    int max_iter = 200;
    char * cluster_num;
    char * iter;
    int k;
    int sign = 0;
    int r;
    /* make vectors linked list */
    struct vector *head_vec;
    struct vector *curr_vec;
    struct cord *head_cord;
    struct cord *curr_cord;
    double n;
    char c;
    int total_cord_num = 0;

    if (argc == 1){
        printf("An Error Has Occurred");
        return 1;
    }

    if (argc > 3){
        printf("An Error Has Occurred");
        return 1;
    }

    if (argc != 2 && argc != 3){
        printf("Invalid Input\n");
        sign++;
    }
    k = strtol(argv[1], &cluster_num, 10);
    if (argc == 3){
        max_iter = strtol(argv[2], &iter, 10);
    }
    if (k <= 1){
        printf("Invalid number of clusters!\n");
        sign++;
    }
    if (max_iter <= 1 || max_iter >= 1000){
        printf("Invalid maximum iteration!\n");
        sign++;
    }
    if(sign > 0){
        free_exit(cluster_num, iter);
        return 1;
    }


    head_cord = malloc(sizeof(struct cord));
    if(head_cord == NULL){
        printf("An Error Has Occurred");
        free_exit(cluster_num, iter);
        return 1;
    }

    curr_cord = head_cord;
    curr_cord->next = NULL;

    head_vec = malloc(sizeof(struct vector));
    if(head_vec == NULL){
        printf("An Error Has Occurred");
        free(head_cord);
        free_exit(cluster_num, iter);
        return 1;
    }
    curr_vec = head_vec;
    curr_vec->next = NULL;


    while (scanf("%lf%c", &n, &c) == 2){
        if (c == '\n'){
            curr_cord->value = n;
            curr_vec->cords = head_cord;
            curr_vec->next = malloc(sizeof(struct vector));
            if(curr_vec->next == NULL){
                printf("An Error Has Occurred");
                free_exit(cluster_num, iter);
                free_LinkedList(head_vec,total_cord_num);
                return 1;
            }
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
            head_cord = malloc(sizeof(struct cord));
            total_cord_num++;
            if(head_cord == NULL){
                printf("An Error Has Occurred");
                free_LinkedList(head_vec,total_cord_num);
                return 1;
            }
            curr_cord = head_cord;
            curr_cord->next = NULL;
            continue;
        }
        curr_cord->value = n;
        curr_cord->next = malloc(sizeof(struct cord));
        if(curr_cord->next == NULL){
            printf("An Error Has Occurred");
            free_LinkedList(head_vec,total_cord_num);
            return 1;
        }
        curr_cord = curr_cord->next;
        curr_cord->next = NULL;
    }


    r = k_means(k, max_iter, head_vec);
    free(head_cord);
    return r;
}
