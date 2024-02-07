import tensorflow as tf
import numpy as np

from keras.models import Sequential
from keras.layers import Dense, Conv1D, Activation
from keras.optimizers import Adam, SGD, RMSprop
from keras.utils import np_utils

model = Sequential()
model.add(Dense(30000, input_dim=62, activation='relu'))
model.add(Dense(8000, activation='relu'))
model.add(Dense(200, activation='relu'))
model.add(Dense(20, activation='relu'))
model.add(Dense(2, activation='softmax'))

opt = RMSprop(lr=0.00000000000000000000000000000000000001)

model.compile(loss='binary_crossentropy',
              optimizer=opt,
              metrics=['accuracy'])

folder = ["Enriched","Enhanced","Elevated"]
file_name = folder[2]

cutoff = 100
min_size = 200

expression = {}

with open(file_name + "TE.csv", "r") as file:
    data = file.readlines()
    for line in data:
        entry = line.strip().split(',')
        tissue = int(entry[0])
        cdb = []
        for i in entry[1:]:
            cdb.append(float(i)*10)
        if tissue in expression.keys():
            expression[tissue].append(np.asarray(cdb))
        else:
            expression[tissue] = []
            expression[tissue].append(np.asarray(cdb))

finished = []
keys = []
models = []
for tissue in expression.keys():
    finished.append(tissue)
    if len(expression[tissue]) >= min_size:
        np.random.shuffle(expression[tissue])
        for tis2 in expression.keys():
            if tis2 not in finished:
                if len(expression[tis2]) >= min_size:
                    train_x = []
                    train_y = []
                    test_x = []
                    test_y = []
                    np.random.shuffle(expression[tis2])
                    for i in range(cutoff):
                        train_x.append(expression[tissue][i])
                        train_x.append(expression[tis2][i])
                        train_y.append(0)
                        train_y.append(1)
                    for i in range(min_size-cutoff):
                        test_x.append(expression[tissue][cutoff + i])
                        test_x.append(expression[tis2][cutoff + i])
                        test_y.append(0)
                        test_y.append(1)                        
                    train_y = np.asarray(train_y)
                    train_y = np_utils.to_categorical(train_y, 2)
                    train_x = np.asarray(train_x)
                    test_y = np.asarray(test_y)
                    test_y = np_utils.to_categorical(test_y, 2)
                    test_x = np.asarray(test_x)
                    keys.append((tissue, tis2))
                    models.append((model.fit(train_x, train_y, batch_size = 200), model.evaluate(test_x, test_y)))
        break     

                
        
