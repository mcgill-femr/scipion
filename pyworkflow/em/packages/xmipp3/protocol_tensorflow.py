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

from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, Activation, InputLayer
from keras.layers import Input, Convolution2D, MaxPooling2D, UpSampling2D
from keras.models import Model
from keras.optimizers import SGD, Adam
from keras.utils import np_utils
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

batch_size = 32
num_classes = 10
epochs = 100

model = Sequential()
model.add(Dense(32, input_dim=100))
model.add(Dense(10,activation='softmax'))
model.add(Activation('relu'))
model.compile(optimizer='rmsprop', loss='categorical_crossentropy', metrics=['accuracy'])

data = np.random.random((1000,100))
labels = np.random.randint(10, size=(1000,1))

one_hot_labels = np_utils.to_categorical(labels,10)
history = model.fit(data, one_hot_labels, epochs=200, batch_size=49)

plt.figure()
plt.plot(history.history['loss'])
plt.show()