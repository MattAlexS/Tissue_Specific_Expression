import tensorflow as tf
import numpy as np

from keras.models import Sequential
from keras.layers import Dense, Conv1D, Activation

folder = ["Enriched","Enhanced","Elevated"]
file_name = folder[1]

factor = 4

expression = {}

with open(file_name + "TE.csv", "r") as file:
    data = file.readlines()
    for line in data:
        entry = line.strip().split(',')
        tissue = int(entry[0])
        cdb = []
        for i in entry[1:]:
            cdb.append(float(i))
        if tissue in expression.keys():
            expression[tissue].append(np.asarray(cdb))
        else:
            expression[tissue] = []
            expression[tissue].append(np.asarray(cdb))

train_y = []
train_x = []
test_y = []
test_x = []

newTis = {}
count = 0

for tissue in expression.keys():
    if len(expression[tissue]) > 199:
        newTis[count] = tissue
        np.random.shuffle(expression[tissue])
        for i in expression[tissue][:180]:
            train_y.append(count)
            train_x.append(i)
        for i in expression[tissue][180:200]:
            test_y.append(count)
            test_x.append(i)
        count += 1

train_y = np.asarray(train_y)
train_x = np.asarray(train_x)
test_y = np.asarray(test_y)
test_x = np.asarray(test_x)

model = Sequential()
model.add(Dense(62, input_dim=62))
model.add(Dense(55, activation='relu'))
model.add(Dense(len(newTis.keys()), activation='softmax'))



model.compile(loss='cosine_proximity',
              optimizer='adam',
              metrics=['accuracy'])

model.fit(train_x, train_y)




        

