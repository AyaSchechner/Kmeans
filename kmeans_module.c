# define PY_SSIZE_T_CLEAN
# include <Python.h>

# include <stdio.h>
# include <stdlib.h>
# include <math.h>


double ** create_vectors_mat(int, int);
double ** create_centroids_mat(int, int);
double euclidean_dis(double *, double *, int);
int check_centroids(double *, int, double);
void calc_new_cent(double **, int *, double **, double *, int, int *, int, int);
void adjust_centroids(double **, double **, int, int *, int *, int, int, double *);
void assign_vectors(double **, double **, int *, int *, int, int, int);
void first_assign_vectors(double **,  double **, int *, int *, int, int, int);
void free_mat(double **, int);
void free_mem(int *, int *, double *);
double ** k_means(int, int, double, double **, double **, int, int);
// static PyObject* fit(PyObject *, PyObject *);


/* creating matrix to store the vectors */
double ** create_vectors_mat(int vec_num, int vec_len){
    double ** vectors_mat;
    int i = 0;
    vectors_mat = (double **) malloc(vec_num * sizeof(double *));
    if(vectors_mat == NULL) {
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
    return vectors_mat;
}


/* creating matrix to store centroids */
double ** create_centroids_mat(int k, int vec_len) {
    double ** centroids_mat;
    int i = 0;
    centroids_mat = (double **) malloc(k * sizeof(double *));
    if (centroids_mat == NULL) {
        printf("An Error Has Occurred");
        return NULL;
    }

    while (i < k) {
        centroids_mat[i] = (double *) malloc(vec_len * sizeof(double));
        if (centroids_mat[i] == NULL) {
            printf("An Error Has Occurred");
            free_mat(centroids_mat, i);
            return NULL;
        }
        i++;
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
int check_centroids(double *centroids_delta, int k, double epsilon){
    int i;
    for (i = 0; i < k; i++) {
        if(centroids_delta[i] >= epsilon){
            return 0;
        }
    }
    return 1;
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


/* assigning vectors to clusters */
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


/* freeing memory */
void free_mat(double ** matrix, int rows_num) {
    int i;
    for (i = 0; i < rows_num; i++) {
        free(matrix[i]);
    }
    free(matrix);
}


/* freeing memory */
void free_mem(int * index, int * cluster_size, double * centroids_deltas){
    free(cluster_size);
    free(centroids_deltas);
    free(index);
}


/* k_means function */
double ** k_means(int k, int iter, double epsilon, double ** vectors_mat, double ** centroids_mat, int vec_num, int vec_len){
    int * index;
    int * cluster_size;
    int check;
    int iter_num;
    double * centroids_deltas;
    if(k >= vec_num){
        printf("Invalid number of clusters!");
        return NULL;
    }

    index = malloc(vec_num * sizeof(int));      /* index[i] = j - vector i belongs to cluster j (size N) */
    if(index == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }

    cluster_size = calloc(k , sizeof(int));    /* cluster_size[i] = j - there are j vectors in cluster i (size k) */
    if(cluster_size == NULL){
        printf("An Error Has Occurred");
        free(index);
        return NULL;
    }

    centroids_deltas = malloc(k * sizeof(double));     /* delta(prev cluster, current cluster) */
    if(centroids_deltas == NULL){
        printf("An Error Has Occurred");
        free(index);
        free(cluster_size);
        return NULL;
    }
    first_assign_vectors(vectors_mat, centroids_mat, index, cluster_size, vec_num, vec_len, k);
    adjust_centroids(vectors_mat, centroids_mat, k, index, cluster_size, vec_num, vec_len, centroids_deltas);
    check = check_centroids(centroids_deltas, k, epsilon);
    if (check == 1){        /* if all deltas are smaller then epsilon */
        free_mem(index, cluster_size, centroids_deltas);
        return centroids_mat;
    }
    iter_num = 1;
    while (iter_num < iter){
        assign_vectors(vectors_mat, centroids_mat, index, cluster_size, vec_num, vec_len, k);
        adjust_centroids(vectors_mat, centroids_mat, k, index, cluster_size, vec_num, vec_len, centroids_deltas);
        check = check_centroids(centroids_deltas, k, epsilon);
        if (check == 1){
            free_mem(index, cluster_size, centroids_deltas);
            return centroids_mat;
        }
        iter_num++;
    }

    free_mem(index, cluster_size, centroids_deltas);
    return centroids_mat;
}


/* kmeans wrap function */
static PyObject* fit(PyObject *self, PyObject *args){
    int k;
    int iter;
    int j;
    int p;
    int i;
    double epsilon;
    int vec_num;
    int vec_len;
    PyObject *lst1;
    PyObject *lst2;
    PyObject *curr_vec;
    PyObject *curr_cent;
    PyObject *answer;
    double** vectors_mat;
    double** centroids_mat;

    if(!PyArg_ParseTuple(args, "iifOO", &k, &iter, &epsilon,  &lst1, &lst2)) {  /* transfer python variables to C*/
        return NULL;
    }

    /* extract vectors number and size */
    vec_num = PyObject_Length(lst1);    
    vec_len = PyObject_Length(PyList_GetItem(lst1, 0));

    /* create vectors matrix and centroids matrix */
    vectors_mat = create_vectors_mat(vec_num, vec_len);
    if(vectors_mat == NULL) {
        return NULL;
    }

    centroids_mat = create_centroids_mat(k, vec_len);
    if(centroids_mat == NULL) {
        free_mat(vectors_mat, vec_num);
        return NULL;
    }

    for (j = 0; j < vec_num; ++j) {
        curr_vec = PyList_GetItem(lst1, j);
        for (p = 0; p < vec_len; ++p){
            vectors_mat[j][p] = PyFloat_AsDouble(PyList_GetItem(curr_vec, p));
        }
    }

    for (j = 0; j < k; ++j) {
        curr_vec = PyList_GetItem(lst2, j);
        for (p = 0; p < vec_len; ++p) {
            centroids_mat[j][p] = PyFloat_AsDouble(PyList_GetItem(curr_vec, p));
        }
    }

    centroids_mat = k_means(k, iter, epsilon, vectors_mat, centroids_mat, vec_num, vec_len);

    /* transfer C final values to Python values*/
    answer = PyList_New(k);
    for (i = 0; i < k; ++i){
        curr_cent = PyList_New(vec_len);
        for (j = 0; j < vec_len; ++j) {
            PyList_SetItem(curr_cent, j, Py_BuildValue("d", centroids_mat[i][j]));
        }
        PyList_SetItem(answer, i, curr_cent);
    }

    free_mat(centroids_mat, k);
    free_mat(vectors_mat, vec_num);
    return answer;
}

static PyMethodDef kmeansMethods[] = {{"fit", (PyCFunction) fit, METH_VARARGS, PyDoc_STR("function expects: number of clusters, max iteration, epsilon, vectors list, centroitds list")}, {NULL,NULL,0,NULL}};

static struct PyModuleDef kmeansmodule = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp", /* name of module */
        NULL, /* module documentation, may be NULL */
        -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        kmeansMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&kmeansmodule);
    if (!m) {
        Py_RETURN_NONE;
    }
    return m;
}
