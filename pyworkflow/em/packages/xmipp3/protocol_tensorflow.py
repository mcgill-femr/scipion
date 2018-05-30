# **************************************************************************
# *
# * Authors:     Javier Mota Garcia (jmota@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import numpy as np
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, BatchNormalization, Activation
from keras.layers import Conv2D, MaxPooling2D, ZeroPadding2D
from keras.optimizers import SGD
from keras.datasets import cifar10
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

num_classes = 10
# Generate dummy data
'''x_train = np.random.random((100, 100, 100, 3))
y_train = keras.utils.to_categorical(np.random.randint(10, size=(100, 1)), num_classes=10)
x_test = np.random.random((20, 100, 100, 3))
y_test = keras.utils.to_categorical(np.random.randint(10, size=(20, 1)), num_classes=10)'''

(x_train, y_train), (x_test, y_test) = cifar10.load_data()

y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)
model = Sequential()
# input: 100x100 images with 3 channels -> (100, 100, 3) tensors.
# this applies 32 convolution filters of size 3x3 each.
model.add(Conv2D(32, (5, 5),padding='same',input_shape=x_train.shape[1:]))
model.add(BatchNormalization())
#model.add(MaxPooling2D(pool_size=(3, 3)))
model.add(Dropout(0.1))
model.add(ZeroPadding2D(padding=(5,5)))
model.add(Conv2D(32, (5, 5)))
model.add(BatchNormalization())
model.add(Activation('leakyrelu'))
model.add(MaxPooling2D(pool_size=(3, 3)))
model.add(Dropout(0.25))

model.add(ZeroPadding2D(padding=(3,3)))
model.add(Conv2D(64, (3, 3)))
model.add(BatchNormalization())
model.add(Dropout(0.1))
#model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(ZeroPadding2D(padding=(3,3)))
model.add(Conv2D(64, (3, 3)))
model.add(BatchNormalization())
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Flatten())
model.add(Dense(256, activation='relu'))
model.add(BatchNormalization())
model.add(Dropout(0.5))
model.add(Dense(10, activation='softmax'))

#sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

history = model.fit(x_train, y_train, validation_split=0.2, batch_size=32, epochs=50)
print history.history.keys()
score = model.evaluate(x_test, y_test, batch_size=32)
print score
'''prediction = model.predict_classes(x_test, batch_size=32)
ok = 0
for i, x in enumerate(prediction):
    if x == np.where(y_test[i]==np.max(y_test[i]))[0]:
        ok += 1

accuracy = ok/float(len(y_test))
print accuracy'''
plt.figure()
plt.plot(history.history['loss'])
plt.figure()
plt.plot(history.history['acc'])
plt.figure()
plt.plot(history.history['val_acc'])
plt.show()
