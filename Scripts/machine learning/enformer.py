# load the data
import numpy as np
import tensorflow as tf
import tensorflow_addons as tfa
import pandas
import joblib

trainingSet=pandas.read_csv('trainingSet_noncoding_v1_moreAnnotations.txt',sep='\t',low_memory=False)
validationSet=pandas.read_csv('validationSet_noncoding_v1_moreAnnotations.txt',sep='\t',low_memory=False)
trainingSet=trainingSet.fillna(0)
validationSet=validationSet.fillna(0)
trainingSet['varID']=trainingSet.loc[:,['hg19_chr','hg19_pos(1-based)','ref','alt']].astype(str).agg('_'.join,axis=1)
validationSet['varID']=validationSet.loc[:,['hg19_chr','hg19_pos(1-based)','ref','alt']].astype(str).agg('_'.join,axis=1)

# get rid of any rows with N's
trainingSet=trainingSet.loc[~((trainingSet['WTSeq'].str.contains('N')) | (trainingSet['AltSeq'].str.contains('N'))),:].reset_index(drop=True)
validationSet=validationSet.loc[~((validationSet['WTSeq'].str.contains('N')) | (validationSet['AltSeq'].str.contains('N'))),:].reset_index(drop=True)


# set up data generator
class DataGenerator(tf.keras.utils.Sequence):
    def __init__(self, list_IDs, labels, dataFrameIn, batch_size=8, top_k_features=32, n_classes=3, shuffle=True):
        self.batch_size = batch_size
        self.labels = labels
        self.list_IDs = list_IDs
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.dataFrameIn=dataFrameIn
        self.top_k_features=top_k_features
        self.on_epoch_end()

    def __len__(self):
        'Denotes the number of batches per epoch'
        if (len(self.list_IDs) % self.batch_size) == 0:
            return int(np.floor(len(self.list_IDs) / self.batch_size))
        else:
            return int(np.ceil(len(self.list_IDs) / self.batch_size))

    def __getitem__(self, index):
        'Generate one batch of data'
        # Generate indexes of the batch
        if (((len(self.list_IDs) % self.batch_size) != 0) & (((index+1)*self.batch_size)>len(self.list_IDs))):
            indexes = self.indexes[index*self.batch_size:]
        else:
            indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]

        # Find list of IDs
        list_IDs_temp = [self.list_IDs[k] for k in indexes]

        # Generate data
        X, y = self.__data_generation(list_IDs_temp)

        return X, y

    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes = np.arange(len(self.list_IDs))
        if self.shuffle == True:
            np.random.shuffle(self.indexes)

    def __data_generation(self, list_IDs_temp):
        'Generates data containing batch_size samples'
        thisBatchSize=len(list_IDs_temp)
        seqs=np.zeros((thisBatchSize,896,self.top_k_features*2))
        for i,ID in enumerate(list_IDs_temp):
            wt=np.load("EnformerEmbeddingsPCA/" + str(self.dataFrameIn.loc[ID,'varID']) + "_WT.npz",allow_pickle=True)
            wt=wt['arr_0']
            alt=np.load("EnformerEmbeddingsPCA/" + str(self.dataFrameIn.loc[ID,'varID']) + "_Alt.npz",allow_pickle=True)
            alt=alt['arr_0']
            # select top K features
            wt=wt[:,:self.top_k_features]
            alt=alt[:,:self.top_k_features]
            seqs[i,:,:self.top_k_features]=wt
            seqs[i,:,self.top_k_features:]=alt
        # non-seq info
        nonSeq = self.dataFrameIn.loc[list_IDs_temp,['pLI','pNull','pRec','oe_mis','oe_lof','oe_mis_upper','oe_lof_upper','AF','nhomalt','GERPN','GERPS','phyloP']].to_numpy()
        #nonSeq = self.dataFrameIn.loc[list_IDs_temp,['pLI','pNull','pRec','oe_mis','oe_lof','oe_mis_upper','oe_lof_upper','GERPN','GERPS']].to_numpy()
        #nonSeq = self.dataFrameIn.loc[list_IDs_temp,['pLI','pNull','pRec','oe_mis','oe_lof','oe_mis_upper','oe_lof_upper']].to_numpy()
        # Store class
        y = self.labels[list_IDs_temp]
        X = {'seqs':seqs,'nonSeqInfo':nonSeq}
        return X, tf.keras.utils.to_categorical(y, num_classes=self.n_classes)


import tensorflow as tf
import tensorflow_addons as tfa
import official.nlp
import official.nlp.keras_nlp.layers

# Construct the classifier as a keras model
def noncodingPathogenicityModel(input_shape=896,numFeatures=32, classes=3,classifier_activation='softmax',**kwargs):
    input0 = tf.keras.layers.Input(shape=(input_shape,numFeatures*2),name='seqs',dtype='float32')
    input1 = tf.keras.layers.Input(shape=12,name='nonSeqInfo')

    seqs = input0
    structured = input1

    posEmb = official.nlp.keras_nlp.layers.PositionEmbedding(max_length=input_shape)
    posEmb_seqs = posEmb(seqs)
    seqs = tf.keras.layers.Add()([seqs,posEmb_seqs])
    seqs = tf.keras.layers.LayerNormalization(axis=-1, epsilon=1e-12, dtype=tf.float32)(seqs)
    
    #seqs = tf.keras.layers.Dense(64)(seqs)

    seqs = official.nlp.keras_nlp.layers.TransformerEncoderBlock(8,256,tf.keras.activations.gelu,output_dropout=0.1,inner_dropout=0.1)(seqs)
    seqs = official.nlp.keras_nlp.layers.TransformerEncoderBlock(8,256,tf.keras.activations.gelu,output_dropout=0.1,inner_dropout=0.1)(seqs)
    seqs = official.nlp.keras_nlp.layers.TransformerEncoderBlock(8,256,tf.keras.activations.gelu,output_dropout=0.1,inner_dropout=0.1)(seqs)
    seqs = official.nlp.keras_nlp.layers.TransformerEncoderBlock(8,256,tf.keras.activations.gelu,output_dropout=0.1,inner_dropout=0.1)(seqs)

    seqs = tf.keras.layers.GlobalMaxPooling1D()(seqs)
    #seqs = tf.keras.layers.Flatten()(seqs)

    x = tf.keras.layers.Concatenate(axis=1)([seqs,structured])
    x = tf.keras.layers.Dense(512,activation='relu')(x)
    x = tf.keras.layers.Dropout(0.1)(x)
    x = tf.keras.layers.Dense(64,activation='relu')(x)
    x = tf.keras.layers.Dropout(0.1)(x)
    x = tf.keras.layers.Dense(8,activation='relu')(x)
    x = tf.keras.layers.Dropout(0.1)(x)
    x = tf.keras.layers.Dense(classes,activation=classifier_activation)(x)

    model = tf.keras.Model(inputs=[input0,input1],outputs=x)

    optimizer = tfa.optimizers.LAMB(learning_rate=1e-3)
    model.compile(optimizer=optimizer,loss='categorical_crossentropy',metrics=['accuracy'])
    return model


batchSize=64

training_generator=DataGenerator(np.arange(len(trainingSet)),trainingSet.loc[:,'classLabel'].values,dataFrameIn=trainingSet,batch_size=batchSize,shuffle=True)
validation_generator=DataGenerator(np.arange(len(validationSet)),validationSet.loc[:,'classLabel'].values,dataFrameIn=validationSet,batch_size=batchSize,shuffle=False)


# model
modelNum=1
numEpochs=5
modelWeightsName='weights_noncodingPathogenicity_model_' + str(modelNum)
callbacks=[tf.keras.callbacks.ModelCheckpoint(filepath=modelWeightsName,save_weights_only=True,save_best_only=True,monitor='val_loss',verbose=1)]

model = noncodingPathogenicityModel()
model.summary()
#model.load_weights(modelWeightsName)
model.fit(training_generator,epochs=numEpochs,callbacks=callbacks,validation_data=validation_generator,class_weight={0:1,1:10,2:5})
model.load_weights(modelWeightsName)

model.evaluate(validation_generator)
y_valid_pred=model.predict(validation_generator)
from sklearn.metrics import classification_report
print("Validation set performance")
print(classification_report(validationSet.loc[:,'classLabel'].to_numpy(), np.argmax(y_valid_pred,axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))

known=pandas.read_csv('knownGenes_noncoding_v1_moreAnnotations.txt',sep='\t',low_memory=False)
novel=pandas.read_csv('novelGenes_noncoding_v1_moreAnnotations.txt',sep='\t',low_memory=False)
known=known.fillna(0)
novel=novel.fillna(0)
known['varID']=known.loc[:,['hg19_chr','hg19_pos(1-based)','ref','alt']].astype(str).agg('_'.join,axis=1)
novel['varID']=novel.loc[:,['hg19_chr','hg19_pos(1-based)','ref','alt']].astype(str).agg('_'.join,axis=1)

# get rid of any rows with N's
known=known.loc[~((known['WTSeq'].str.contains('N')) | (trainingSet['AltSeq'].str.contains('N'))),:].reset_index(drop=True)
novel=novel.loc[~((novel['WTSeq'].str.contains('N')) | (novel['AltSeq'].str.contains('N'))),:].reset_index(drop=True)

known_generator=DataGenerator(np.arange(len(known)),known.loc[:,'classLabel'].values,dataFrameIn=known,batch_size=batchSize,shuffle=False)
novel_generator=DataGenerator(np.arange(len(novel)),novel.loc[:,'classLabel'].values,dataFrameIn=novel,batch_size=batchSize,shuffle=False)


model.evaluate(known_generator)
y_known_pred=model.predict(known_generator)
print("Known genes set performance")
print(classification_report(known.loc[:,'classLabel'].to_numpy(), np.argmax(y_known_pred,axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))
print("Know genes splicing")
print(classification_report(known.loc[known['varType']=='splicing','classLabel'].to_numpy(), np.argmax(y_known_pred[known['varType']=='splicing'],axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))
print('Know genes intronic')
print(classification_report(known.loc[known['varType']=='intronic','classLabel'].to_numpy(), np.argmax(y_known_pred[known['varType']=='intronic'],axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))
print('Know genes UTR3')
print(classification_report(known.loc[known['varType']=='UTR3','classLabel'].to_numpy(), np.argmax(y_known_pred[known['varType']=='UTR3'],axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))
print('Know genes UTR5')
print(classification_report(known.loc[known['varType']=='UTR5','classLabel'].to_numpy(), np.argmax(y_known_pred[known['varType']=='UTR5'],axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))
print('Know genes upstream')
print(classification_report(known.loc[known['varType']=='upstream','classLabel'].to_numpy(), np.argmax(y_known_pred[known['varType']=='upstream'],axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))

known['classLabelBinary']=0
known.loc[known['classLabel']>0,'classLabelBinary']=1
print("Known genes set performance (binary)")
print(classification_report(known.loc[:,'classLabelBinary'].to_numpy(), np.round(1-y_known_pred[:,0]).astype('int'), target_names=['Benign','Pathogenic'], digits=3))



model.evaluate(novel_generator)
y_novel_pred=model.predict(novel_generator)
print("Novel genes set performance")
print(classification_report(novel.loc[:,'classLabel'].to_numpy(), np.argmax(y_novel_pred,axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))

novel['classLabelBinary']=0
novel.loc[novel['classLabel']>0,'classLabelBinary']=1
print("Novel genes set performance (binary)")
print(classification_report(novel.loc[:,'classLabelBinary'].to_numpy(), np.round(1-y_novel_pred[:,0]).astype('int'), target_names=['Benign','Pathogenic'], digits=3))

training_generator=DataGenerator(np.arange(len(trainingSet)),trainingSet.loc[:,'classLabel'].values,dataFrameIn=trainingSet,batch_size=batchSize,shuffle=False)
model.evaluate(training_generator)
training_generator
y_training_pred=model.predict(training_generator)
print("training genes set performance")
print(classification_report(training.loc[:,'classLabel'].to_numpy(), np.argmax(y_training_pred,axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))
print("training genes splicing")
print(classification_report(trainingSet.loc[trainingSet['varType']=='splicing','classLabel'].to_numpy(), np.argmax(y_training_pred[trainingSet['varType']=='splicing'],axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))
print('training genes intronic')
print(classification_report(trainingSet.loc[trainingSet['varType']=='intronic','classLabel'].to_numpy(), np.argmax(y_training_pred[trainingSet['varType']=='intronic'],axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))
print('training genes UTR3')
print(classification_report(trainingSet.loc[trainingSet['varType']=='UTR3','classLabel'].to_numpy(), np.argmax(y_training_pred[trainingSet['varType']=='UTR3'],axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))
print('training genes UTR5')
print(classification_report(trainingSet.loc[trainingSet['varType']=='UTR5','classLabel'].to_numpy(), np.argmax(y_training_pred[trainingSet['varType']=='UTR5'],axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))
print('training genes upstream')
print(classification_report(trainingSet.loc[trainingSet['varType']=='upstream','classLabel'].to_numpy(), np.argmax(y_training_pred[trainingSet['varType']=='upstream'],axis=1).astype('int'), target_names=['Benign','Dominant','Recessive'], digits=3))

