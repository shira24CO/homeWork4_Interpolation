''' Ilan Shklover 206753550
Shira Cohen 211777834
Eli Ben Aharon 311614853
'''

All_Elementary_matrix= {}

def liniar(p1,p2, x,point):
    f = ((point[p1]-point[p2])/(p1-p2))*x+(p1*point[p2]-p2*point[p1])/(p1-p2)
    return f
def polinomit(point,x):
    matrix =[[0 for j in range(len(point))]for i in range(len(point))]
    keys = list(point.keys())
    for j in range(len(point)):
        for i in range(len(point)):
            if j == 0:
                matrix[i][0] = 1
            else:
                matrix[i][j] = pow(keys[i],j)
    b = [[0 for j in range(1)]for i in range(len(point))]
    for i in range(len(matrix)):
        b[i][0] = point[keys[i]]
    result = elementary_matrix(matrix,b)
    f = 0
    for k in range(len(result)):
        f += result[k][0]*pow(x,k)
    return f

def lagrange(point, x):
    keys = list(point.keys())
    L_I = 1
    sum = 0
    for i in range(len(point)):
        for j in range(len(point)):
            if i != j:
                L_I *= ((x-keys[j]))/(keys[i]-keys[j])
        sum += L_I * point[keys[i]]
        L_I = 1
    return sum


def nevil(point,x):
    keys = list(point.keys())
    keyspnm = list(point.keys())
    y= []
    for i in range(len(keyspnm)):
        y.append(point[ keyspnm[i]])
    pnm = y.copy()
    counter = 1
    for i in range(len(keyspnm) - 1):
        temp =[]
        for j in range(len(keyspnm)-counter):
            temp.append(((x- keyspnm[j])*pnm[j+1]-(x -keyspnm[j+counter])*pnm[j])/(keyspnm[j+counter]-keyspnm[j]))
        pnm = temp.copy()
        counter+=1
    return pnm[0]


def spline(point, x):
    keys = list(point.keys())
    y = [point[keys[i]] for i in range(len(point))]
    h = [(keys[i + 1] - keys[i]) for i in range(len(point)-1)]
    l = [h[i]/(h[i-1]+h[i]) for i in range(1,len(h))]
    m = [1-l[i] for i in range(len(l))]
    d = [[(6/(h[i-1]+h[i]))*(((y[i+1]-y[i])/h[i])-((y[i]-y[i-1])/h[i-1]))] for i in range(1,len(h))]
    matrix = create_I_matrix(len(point)-2)
    for i in range(len(matrix)): # הצבת 2
        matrix[i][i] = 2
    for j in range(1,len(matrix)) : # הצבת ה m
        for k in range(len(matrix)-1):
            matrix[j][k] = m[k+1]
    for line in range(len(matrix)-1):  # הצבת ה l
        for col in range(1,len(matrix)):
            matrix[line][col] = l[line]
    M = elementary_matrix(matrix,d)
    lower_p = 0
    for t in range(len(keys)):
        if keys[t] < x:
            lower_p = t
    M = [[0]]+M+[[0]]  # the second derivative of the points in the upper and lower bounderies are equal to 0.
    a = ((((pow((keys[lower_p+1]-x),3))*M[lower_p][0])+((pow((x-keys[lower_p]),3))*M[lower_p+1][0]))/(6*h[lower_p]))
    b = (((keys[lower_p+1]-x)*y[lower_p])+((x-keys[lower_p])*y[lower_p+1]))/h[lower_p]
    c = ((((keys[lower_p+1]-x)*M[lower_p][0])+((x-keys[lower_p])*M[lower_p+1][0]))*h[lower_p])/6
    s = a+b-c
    return s




def printmat(matrix):
    print("\n")
    for i in matrix:
        for j in i:
            print(f'{j}  ')
        print('\n')
    print('\n')






def product_calculation(vector_line, vector_column):  # Multiply line and column
    result = 0
    for element in range(len(vector_line)):
        result += vector_line[element] * vector_column[element]
    return result


def swap_lines_of_matrix(matrix, index_line1, index_line2):
    for column in range(len(matrix)):
        temp_value = matrix[index_line2][column]
        matrix[index_line2][column] = matrix[index_line1][column]
        matrix[index_line1][column] = temp_value


def check_column_in_matrix(matrix, some_col):
    for line in range(len(matrix)):
        if matrix[line][some_col] != 0:
            return True  # The column is not full of 0's
    return False  # The column is full of 0's


def swap_columns_of_matrix(matrix, index_column1):  # 'matrix' is a matrix that has a column full of 0's.
    last_column = len(matrix) - 1
    for column in range(last_column, index_column1, 1):
        if check_column_in_matrix(matrix, column) == True:
            for line in range(len(matrix)):
                matrix[line][index_column1] = matrix[line][last_column]
            for index in range(len(matrix)):
                matrix[index][last_column] = 0
        return


def multiply_matrix(mat1, mat2):  # matrixes of n*n
    new_mat = []  # contains the result of the multiplication
    line_index = 0
    vector_col = []
    for line_mat1 in mat1:
        temp_line_mat = []
        for column in range(len(mat2[0])):
            vector_col = []
            for line in range(len(mat1)):
                vector_col.append(mat2[line][column])
            temp_line_mat.append(product_calculation(line_mat1, vector_col))
        new_mat.append(temp_line_mat)
    return new_mat


def create_I_matrix(size):
    matrixI = []
    for i in range(size):
        matrixI_helper = []
        for j in range(size):
            if i == j:
                matrixI_helper.append(1)
            else:
                matrixI_helper.append(0)
        matrixI.append(matrixI_helper)
    return matrixI


def find_max_of_column(matrix, j):
    element = create_I_matrix(len(matrix))
    #  Find the maximun value in a column in order to change lines
    maxElem = abs(matrix[j][j])
    maxRow = j
    for k in range(j + 1, len(matrix)):  # Interacting over the next line,in the same column
        if (abs(matrix[k][j]) > maxElem):
            maxElem = abs(matrix[k][j])  # Next line on the diagonal
            maxRow = k
    swap_lines_of_matrix(element, maxRow, j)
    return element


def print_state(elementary, matrix):
    return multiply_matrix(elementary, matrix)


def elementary_matrix(matrix, result_vector):
    counter_for_elementary_matrix = 0
    counter_for_elementary_operations1 = (pow(len(matrix), 2) + len(
        matrix)) / 2  # In order to create an upper triangular form for matrix 3 *3 we operate 3+2+1 operations(sum of arithmetic progression)
    while counter_for_elementary_matrix != counter_for_elementary_operations1:
        for column in range(len(matrix)):  #
            elementary_mat = find_max_of_column(matrix, column)
            matrix = print_state(elementary_mat,matrix)
            result_vector = multiply_matrix(elementary_mat,result_vector)
            for line in range(len(matrix)):
                if line == column and matrix[line][column] != 0:
                    piv = 1 / matrix[line][column]
                    elementary_mat = create_I_matrix(len(matrix))
                    elementary_mat[line][column] = piv
                    result_vector = multiply_matrix(elementary_mat,result_vector)
                    matrix = print_state(elementary_mat, matrix)
                    counter_for_elementary_matrix += 1
                    All_Elementary_matrix[
                        counter_for_elementary_matrix] = elementary_mat  # Enter new elementary matrix in the dictionary.
                elif line == column and matrix[line][column] == 0:  # we need to swap lines
                    line_to_swap_with = -1
                    for l in range(len(matrix)):
                        if matrix[l][column] != 0:
                            line_to_swap_with = l
                            swap_lines_of_matrix(matrix, line_to_swap_with, line)
                    if line_to_swap_with == -1:  # we did not find in the column 'column' a number different than zero. Therefore we can not swap lines. So,we will try to swap columns.
                        swap_columns_of_matrix(matrix, column)
                elif line != column and line > column:
                    elementary_mat = create_I_matrix(len(matrix))
                    piv = - matrix[line][column] / matrix[column][column]
                    elementary_mat[line][column] = piv
                    matrix = print_state(elementary_mat,matrix)
                    result_vector = multiply_matrix(elementary_mat,result_vector)
                    counter_for_elementary_matrix += 1
                    All_Elementary_matrix[counter_for_elementary_matrix] = elementary_mat
    # Until here we receive an upper triangle matrix
    counter_for_elementary_operations2 = ((pow(len(matrix), 2) + len(matrix)) / 2) - len(matrix)
    counter_for_elementary_matrix2 = 0
    while counter_for_elementary_matrix2 != counter_for_elementary_operations2:
        for column in range(len(matrix) - 1, -1, -1):
            for line in range(column - 1, -1, -1):
                if line != column and line < column:
                    elementary_mat = create_I_matrix(len(matrix))
                    piv = - matrix[line][column] / matrix[column][column]
                    elementary_mat[line][column] = piv
                    matrix = print_state(elementary_mat,matrix)
                    result_vector = multiply_matrix(elementary_mat,result_vector)
                    counter_for_elementary_matrix2 += 1
                    All_Elementary_matrix[
                        counter_for_elementary_matrix + counter_for_elementary_matrix2] = elementary_mat
    return result_vector



point= {1:0, 1.2: 0.11246, 1.3: 0.167996,1.4:0.222709}
keys = list(point.keys())
x = 1
print("liniar", liniar(keys[0],keys[3],x,point))
print("polinomit", polinomit(point,x))
print("lagrange", lagrange(point,x))
print("nevil", nevil(point,x))
print("spline", spline(point,x))