import pandas as pd
import numpy as np
import sys
import math
import mykmeanssp


# create a sorted vectors table from the input files via inner join
def create_table(input_file1, input_file2):
    df1 = pd.read_csv(input_file1, header=None)
    df2 = pd.read_csv(input_file2, header=None)
    df_join = pd.merge(df1, df2, left_on=0, right_on=0, how='inner')
    df_sorted = df_join.sort_values(by=0)
    return df_sorted


# computes Euclidean distance between two vectors
def euclidean_dis(vec1, vec2):
    d = len(vec1)
    dis = 0
    for i in range(d):
        dis += ((vec1[i] - vec2[i]) ** 2)
    dis = (math.sqrt(dis))
    return dis


# calculate distance between a given vector to all centroids
def calculate_dist(vectors, centroids, vecs_dist, cent_index):
    for i in range(len(vectors)):
        check = False
        for j in cent_index:
            if i == j:
                check = True
        if check:
            vecs_dist[i] = 0
            continue
        arg1 = vectors[i]
        arg2 = centroids[0]
        first_dist = euclidean_dis(arg1, arg2)
        for j in range(1, len(centroids)):
            arg3 = vectors[i]
            arg4 = centroids[j]
            second_dist = euclidean_dis(arg3, arg4)
            if second_dist < first_dist:
                first_dist = second_dist
        vecs_dist[i] = first_dist


# choose a centroid (for centroids initialization)
def choose_vec(vectors, vecs_dist):
    dist_sum = 0
    for i in vecs_dist:
        dist_sum += i
    vec_index = [i for i in range(len(vectors))]
    vec_prob = [0 for i in range(len(vectors))]
    for i in range(len(vectors)):
        vec_prob[i] = (vecs_dist[i]) / dist_sum
    chosen_index = np.random.choice(vec_index, p=vec_prob)
    return chosen_index


# create the vectors list
def create_mat1(vectors):
    vectors_mat = []
    for i in range(vectors.shape[0]):
        curr_vec = []
        for j in range(1, vectors.shape[1]):
            curr_vec.append(vectors[vectors[0] == float(i)].values[0][j])
        vectors_mat.append(curr_vec)
    return vectors_mat


# create the centroids list
def create_mat2(centroids):
    cent_mat = []
    for i in centroids:
        curr_vec = []
        for j in range(1, centroids[0][0].size):
            curr_vec.append(i[0][j])
        cent_mat.append(curr_vec)
    return cent_mat


# print the initial centroid's index and the final centroids
def print_centroids(initial_indx, centroids, d):
    for i in range(len(initial_indx)):
        if i < (len(initial_indx) - 1):
            print(initial_indx[i], ",", end="")
        else:
            print(initial_indx[i], '\n', end="")
    for vector in centroids:
        centroid_len = d
        for i in range(d):
            if centroid_len > 1:
                print('{:.4f}'.format(vector[i]), ",", end="")
                centroid_len -= 1
            else:
                print('{:.4f}'.format(vector[i]), '\n', end="")


# k_means_pp function
def k_means_pp(k, iter_num, eps, input_file1, input_file2):
    cent_index = []
    centroids_num = 0
    sorted_vec = create_table(input_file1, input_file2)
    vec_num = sorted_vec.shape[0]
    if k >= vec_num:
        print("Invalid number of clusters!")
        exit()
    vec_len = sorted_vec.shape[1] - 1
    np.random.seed(0)
    rand = float(np.random.choice(sorted_vec.shape[0]))
    first_centroid = sorted_vec[sorted_vec[0] == rand].values
    cent_index.append(int(rand))
    centroids_mat = [first_centroid]
    centroids_num += 1
    sorted_vec2 = create_table(input_file1, input_file2)
    cen_mat = create_mat2(centroids_mat)
    vec_mat = create_mat1(sorted_vec2)
    while centroids_num < k:
        vecs_dist = [0 for i in vec_mat]
        calculate_dist(vec_mat, cen_mat, vecs_dist, cent_index)
        chosen_one = choose_vec(vec_mat, vecs_dist)
        cent_index.append(chosen_one)
        curr_vec = vec_mat[chosen_one]
        cen_mat.append(curr_vec)
        sorted_vec = sorted_vec[sorted_vec[0] != chosen_one]
        centroids_num += 1
    final_centroids = mykmeanssp.fit(k, iter_num, eps, vec_mat, cen_mat)
    if final_centroids == None:
        print("An Error Has Occurred")
        exit()
    print_centroids(cent_index, final_centroids, vec_len)
    return 1


# extract argumants + basic check ups 
args = sys.argv
args_len = len(args)
max_iter = 300
offset = 0
k1 = 0
epsilon = 0
file1 = None
file2 = None
if args_len != 6 and args_len != 5:
    print("An Error Has Occurred")
    exit()
try:
    k1 = int(args[1])
    if args_len == 6:
        offset = 1
        max_iter = int(args[2])
    epsilon = float(args[2 + offset])
    file1 = args[3 + offset]
    file2 = args[4 + offset]
    if k1 <= 1:
        print("Invalid number of clusters!")
        exit()
    if max_iter <= 1 or max_iter >= 1000:
        print("Invalid maximum iteration!")
        exit()
except Exception:
    print("An Error Has Occurred")
    exit()

k_means_pp(k1, max_iter, epsilon, file1, file2)
