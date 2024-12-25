import numpy as np

# ==== Functions ====
def make_vector():
    vector = np.array([1, 2, 3, 4, 5])
    print("vector:")
    print(vector, "\n")

def make_matrix():
    matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    print("matrix:")
    print(matrix, "\n")

    identity_matrix = np.eye(3)
    print("identity_matrix:")
    print(identity_matrix, "\n")

    zero_matrix = np.zeros((3, 3))
    print("zero_matrix:")
    print(zero_matrix, "\n")

    one_matrix = np.ones((3, 3))
    print("one_matrix:")
    print(one_matrix, "\n")

    random_matrix = np.random.rand(3, 3)
    print("random_matrix:")
    print(random_matrix, "\n")

def make_tensor():
    tensor = np.array([[[1, 2, 3], [4, 5, 6], [7, 8, 9]], [[10, 11, 12], [13, 14, 15], [16, 17, 18]], [[19, 20, 21], [22, 23, 24], [25, 26, 27]]])
    print("tensor:")
    print(tensor, "\n")

def calculate_matrix():
    m1 = np.array([[1, 2], [3, 4]])
    print("m1:")
    print(m1, "\n")

    print("m1 + 3:")
    print(m1 + 3, "\n")
    print("m1 - 3:")
    print(m1 - 3, "\n")
    print("m1 * 3:")
    print(m1 * 3, "\n")
    print("m1 / 3:")
    print(m1 / 3, "\n")

    m2 = np.array([[4, 3], [2, 1]])
    m3 = m1 + m2

    print("m1:")
    print(m1, "\n")
    print("m2:")
    print(m2, "\n")
    print("m1 + m2:")
    print(m3, "\n")

    m4 = np.dot(m1, m2)
    print("np.dot(m1, m2):")
    print(m4, "\n")

    print("m1 @ m2")
    print(m1 @ m2, "\n")

    print("m1 * m2")
    print(m1 * m2)
    print("各要素の積となる\n")

def transform_matrix():
    m1 = np.array([[1, 2], [3, 4]])
    print("m1")
    print(m1, "\n")

    print("m1.T")
    print(m1.T, "\n")

    print("m1.transpose()")
    print(m1.transpose(), "\n")

    print("m1.reshape(4, 1)")
    print(m1.reshape(4, 1), "\n")

def linear_algebra():
    m1 = np.array([[1, 2], [3, 4]])
    print("m1")
    print(m1, "\n")

    print("np.linalg.inv(m1)")
    print(np.linalg.inv(m1), "\n")

    print("m1 @ np.linalg.inv(m1)")
    print(m1 @ np.linalg.inv(m1), "\n")

    print("m1")
    print(m1, "\n")

    eigenvalues, eigenvectors = np.linalg.eig(m1)
    print("eigenvalues:")
    print(eigenvalues, "\n")
    print("eigenvectors:")
    print(eigenvectors, "\n")

    determinant = np.linalg.det(m1)
    print("determinant of m1:")
    print(determinant, "\n")

# ==== Main ====
def main():
    # make_vector()
    # make_matrix()
    # make_tensor()
    # calculate_matrix()
    # transform_matrix()
    linear_algebra()

# スクリプトとして実行された場合にmain関数を呼び出す
if __name__ == "__main__":
    main()
