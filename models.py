from keras.models import Sequential, Model
from keras.layers import Convolution1D, MaxPooling1D, GlobalMaxPooling1D, Dense, Dropout, Flatten, Input, concatenate, BatchNormalization, Activation, add, Bidirectional, LSTM, GRU
from keras import regularizers

# DeepBind model
def DeepBind(shape = None, params = None, penalty = 0.005):
    model = Sequential()
    model.add(Convolution1D(filters=16, kernel_size=13, padding='same', activation='relu',
                kernel_regularizer=regularizers.l2(penalty),input_shape=shape))
    model.add(GlobalMaxPooling1D())
    model.add(Dense(units=32, activation='relu', kernel_regularizer=regularizers.l2(penalty)))
    model.add(Dropout(params['DROPOUT']))
    model.add(Dense(units=1))
    
    print model.summary()
    return model

# build hybrid model
def SharedDeepBindwithShape(shape1=None, shape2=None, params=None, penalty = 0.005):
    
    digit_input = Input(shape=shape1)
    X = Convolution1D(16, 13, activation='relu', padding='same')(digit_input)
    out = GlobalMaxPooling1D()(X)
    share_model = Model(digit_input, out)
    
    main_input = Input(shape=shape1, name='sequence')
    out_main = share_model(main_input)
    
    auxiliary_input = Input(shape=shape2, name='shape')
    auxiliary_conv1 = Convolution1D(4, 1, activation='relu', padding='same', name='shape_conv')(auxiliary_input) 
    out_aux = share_model(auxiliary_conv1)
    
    concat = concatenate([out_main, out_aux], axis=-1)
    Y = Dense(32, activation='relu', kernel_regularizer=regularizers.l2(penalty))(concat)
    Y = Dropout(params['DROPOUT'])(Y)
    output = Dense(1)(Y)
    
    model = Model(inputs=[main_input, auxiliary_input], outputs=output)
    print model.summary()
    return model

# build hybrid model
def SharedDeepCNNwithShape(shape1=None, shape2=None, params=None, penalty = 0.005):
    
    digit_input = Input(shape=shape1)
    X = Convolution1D(16, 13, activation='relu', padding='same')(digit_input)
    X = MaxPooling1D(2, 2)(X)
    X = Dropout(0.2)(X)
    X = Convolution1D(32, 7, activation='relu', padding='same')(X)
    X = MaxPooling1D(2, 2)(X)
    X = Dropout(0.2)(X)
    X = Convolution1D(32, 5, activation='relu', padding='same')(X)
    out = GlobalMaxPooling1D()(X)
    share_model = Model(digit_input, out)
    
    main_input = Input(shape=shape1, name='sequence')
    out_main = share_model(main_input)
    
    auxiliary_input = Input(shape=shape2, name='shape')
    auxiliary_conv1 = Convolution1D(4, 1, activation='relu', padding='same', name='shape_conv')(auxiliary_input) 
    out_aux = share_model(auxiliary_conv1)
    
    concat = concatenate([out_main, out_aux], axis=-1)
    Y = Dense(32, activation='relu', kernel_regularizer=regularizers.l2(penalty))(concat)
    Y = Dropout(params['DROPOUT'])(Y)
    output = Dense(1)(Y)
    
    model = Model(inputs=[main_input, auxiliary_input], outputs=output)
    print model.summary()
    return model

    
# build other models
def DeepCNN(shape = None, params = None, penalty = 0.005):
    model = Sequential()
    model.add(Convolution1D(16, 13, activation='relu', padding='same', input_shape=shape))
    model.add(MaxPooling1D(2, 2))
    model.add(Dropout(0.2))
    model.add(Convolution1D(32, 7, activation='relu', padding='same'))
    model.add(MaxPooling1D(2, 2))
    model.add(Dropout(0.2))
    model.add(Convolution1D(32, 5, activation='relu', padding='same'))
    model.add(GlobalMaxPooling1D())
    model.add(Dense(32, activation='relu', kernel_regularizer=regularizers.l2(penalty)))
    model.add(Dropout(params['DROPOUT']))
    model.add(Dense(1))
    
    print model.summary()
    return model


