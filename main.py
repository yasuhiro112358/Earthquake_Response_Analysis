import numpy as np

vector = np.array([1, 2, 3, 4, 5])
print("vector:")
print(vector, "\n")

matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
print("matrix:")
print(matrix, "\n")

tensor = np.array([[[1, 2, 3], [4, 5, 6], [7, 8, 9]], [[10, 11, 12], [13, 14, 15], [16, 17, 18]], [[19, 20, 21], [22, 23, 24], [25, 26, 27]]])
print("tensor:")
print(tensor, "\n")

print("matrix + 3:")
print(matrix + 3, "\n")
print("matrix - 3:")
print(matrix - 3, "\n")
print("matrix * 3:")
print(matrix * 3, "\n")
print("matrix / 3:")
print(matrix / 3, "\n")

m1 = np.array([[1, 2], [3, 4]])
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

print("m1")
print(m1, "\n")
print("m1.T")
print(m1.T, "\n")
print("m1.transpose()")
print(m1.transpose(), "\n")
print("m1.reshape(4, 1)")
print(m1.reshape(4, 1), "\n")

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

print(eigenvectors[0,1], "\n")

det = np.linalg.det(m1)
print("det:")
print(det, "\n")

identity_matrix = np.eye(3)
print("identity_matrix:")
print(identity_matrix, "\n")

zero_matrix = np.zeros((3, 3))
print("zero_matrix:")
print(zero_matrix, "\n")