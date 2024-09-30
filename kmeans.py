import math
import sys


def extract_vectors(inpt):
    with open(inpt, "r") as file:
        vectors = []
        for line in file:
            if len(line) > 1:
                vectors.append([float(coordinate) for coordinate in line[:len(line)].split(',')])
        return vectors


def euclidean_dis(vec1, vec2):
    d = len(vec1)
    dis = 0
    for i in range(d):
        a = vec1[i]
        b = vec2[i]
        dis += ((vec1[i]-vec2[i])**2)
    dis = (math.sqrt(dis))
    return dis


def assign_centroids(clust, vec):
    for i in range(k):
        clust.append([vec[i]])


def make_centroids(vec):
    centroids = []
    for i in range(k):
        centroids.append(vec[i])
    return centroids


def assign_vectors(clust, vec, centroids):
    for i in range(len(vec)):
        dis_from_first_cent = euclidean_dis(vec[i], centroids[0])
        closest_cent = 0
        for j in range(1, len(centroids)):
            dis_from_curr_cent = euclidean_dis(vec[i], centroids[j])
            if dis_from_curr_cent < dis_from_first_cent:
                dis_from_first_cent = dis_from_curr_cent
                closest_cent = j
        clust[closest_cent].append(vec[i])


def assign_vectors2(vectors, centroids):
    new_cluster = [[] for i in range(k)]
    for i in vectors:
        dis_from_first_cent = euclidean_dis(i, centroids[0])
        closest_cent = 0
        for j in range(1, len(centroids)):
            dis_from_curr_cent = euclidean_dis(i, centroids[j])
            if dis_from_curr_cent < dis_from_first_cent:
                dis_from_first_cent = dis_from_curr_cent
                closest_cent = j
        new_cluster[closest_cent].append(i)
    return new_cluster


def calc_new_cent(cluster):
    new_cent = [0 for i in range(len(cluster[0]))]
    for i in cluster:
        for j in range(len(cluster[0])):
            new_cent[j] += i[j]
    for p in range(len(new_cent)):
        c = new_cent[p]
        d = len(cluster)
        new_cent[p] = ((new_cent[p]) / len(cluster))
        t = new_cent[p]
        f = new_cent[p]
    return new_cent


def adjust_centroids(clusters, centroids):
    delta = [0 for i in range(k)]
    for i in range(k):
        new_cent = calc_new_cent(clusters[i])
        delta[i] = euclidean_dis(centroids[i], new_cent)
        centroids[i] = new_cent
    return delta


def check_centroids(centroids):
    for i in centroids:
        if i >= 0.001:
            return False
    return True


def print_centroids(centroids, d):
    # change ahead:
    for vector in centroids:
        centroid_len = d
        # print("len= ", centroid_len)
        for i in range(centroid_len):
            if centroid_len > 1:
                print('{:.4f}'.format(vector[i]) + ",", end="")
                centroid_len -= 1
            else:
                print('{:.4f}'.format(vector[i]) + '\n', end="")


def k_means(K, iters, input_data):
    vectors = extract_vectors(input_data)
    iter_num = 0
    k1 = int(K)  # k->k1
    centroids = make_centroids(vectors)
    clusters = [[] for i in range(k1)]
    assign_vectors(clusters, vectors, centroids)
    deltas = adjust_centroids(clusters, centroids)
    check = check_centroids(deltas)
    if check:
        for i in centroids:
            print(i)
        exit()
    while iter_num < iters:
        clusters = assign_vectors2(vectors, centroids)
        deltas = adjust_centroids(clusters, centroids)
        check = check_centroids(deltas)
        if check:
            print_centroids(centroids, len(centroids[0]))
            exit()
        iter_num += 1
    print_centroids(centroids, len(centroids[0]))
    exit()


args = sys.argv
args_len = len(args)
max_iter = 200
offset = 0
if args_len != 4 and args_len != 3:
    print("Invalid Input")
    exit()
try:
    k = int(args[1])
    if args_len == 4:
        offset = 1
        max_iter = int(args[2])
    input_file = args[2 + offset]
    if k <= 1:
        print("Invalid number of clusters!")
        exit()
    if max_iter <= 1 or max_iter >= 1000:
        print("Invalid maximum iteration!")
        exit()
except Exception:
    print("Invalid Input")
    exit()


k_means(k, max_iter, input_file)
