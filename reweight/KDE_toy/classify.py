import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import tensorflow as tf
from tensorflow import keras
from keras import backend as K 
from tensorflow.keras import layers

devices = tf.config.list_physical_devices()
print(devices)
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            print(gpu)
            tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
        print('Error', e)

print('##############################################################')


def dataSelP(data, binEdge, weights = None, binEdgeSec = None, weightVal = 1000):
    sel = np.asarray(np.where((data['m4j']>binEdge[0]) & (data['m4j']<=binEdge[1]), True, False))
    if binEdgeSec is not None:
        sel |= np.asarray(np.where((data['m4j']>binEdgeSec[0]) & (data['m4j']<=binEdgeSec[1]), True, False))
    if weights is not None:
        sel &= np.asarray(np.where((weights>weightVal), True, False))
    return sel

seed = 3
sigma = 30
# m4jBinEdges = np.array([[0, 361], [361, 425], [425, 479], [479, 533], [533, 591], [591, 658], [658, 741], [741, 854], [854, 1044], [1044, 1800]]);

m4jBinEdges = np.array([[0, 271], [271, 359], [359, 423], [423, 478], [478, 530], [530, 584], [584, 644], [644, 718], [718, 831], [831, 1800]])
# nBins = 100
data3b = np.load('data3b_toy_seed_'+str(seed)+'.npy', allow_pickle=True).item()
nBins = 15
w3to4FileName = f'hist_sig_w3to4_toy_seed_3_bins_{nBins}.npy'
w3to4 = np.load(w3to4FileName)

data4bSig = np.load('data4b_sig_toy_seed_'+str(seed)+'.npy', allow_pickle=True).item()
w3to4FileNameSig = f'hist_sig_w3to4_toy_seed_3_bins_{nBins}.npy'
w3to4Sig = np.load(w3to4FileNameSig)

data3b['passHLT'], data4bSig['passHLT'] = np.ones(len(data3b['m4j'])), np.ones(len(data4bSig['m4j']))

w3to4Sig[np.where(np.isnan(w3to4Sig) == True)[0]] = 0
w3to4Sig[np.where(np.isinf(w3to4Sig) == True)[0]] = 0

print('Data obtained')

def get_rand_pair(data, binInd = None, weight = None):
    if binInd is None:
        SRsel = np.full(len(data['m4j']), True)
    else:    
        SRsel = dataSelP(data, m4jBinEdges[binInd], weightVal = 0)
    randind = np.random.choice(len(data['m4j'][SRsel]), size=int(len(data['m4j'][SRsel])/2), replace=False)
    data1 = np.array([data['leadStM'][SRsel][randind],  data['sublStM'][SRsel][randind],  data['m4j'][SRsel][randind] ]).T
    data2 = np.array([data['leadStM'][SRsel][~randind], data['sublStM'][SRsel][~randind], data['m4j'][SRsel][~randind]]).T
    if weight is not None:
        w1, w2 = weight[SRsel][randind], weight[SRsel][~randind]
    else:
        w1, w2 = np.ones(len(data1)), np.ones(len(data2))
    return data1, data2, w1, w2

def make_class_plot(pp, b_pb, d_pb, wp0, wp1, binInd = 0, title = '', xlab = '' ):
    figf = plt.figure(figsize = (8,5))
    lo = np.min([np.min(b_pb), np.min(d_pb)])
    hi = np.max([np.max(b_pb), np.max(d_pb)])

    # lo = np.min([0, 1])
    # hi = np.max([0, 1])
    a = plt.hist(b_pb, range = (lo, hi), bins = 30, weights = wp0, label = 'bkg model')
    b = plt.hist(d_pb ,fill = False, hatch = '/', range = (lo, hi), bins = 30,  weights = wp1, label = 'Data');
    plt.legend(); plt.title(f'({title}) 3D classification: SR {binInd}'); 
    plt.xlabel(xlab);
    plt.yscale('log')
    plt.ylabel('counts'); 
    # plt.show(); 
    pp.savefig(figf)
    plt.close()
    return pp

def classify(binInd, do_train = False, do_stat = False):
    np.random.seed(seed)
    tf.random.set_seed(seed)
    tf.keras.utils.set_random_seed(seed)
    # binInd = 5
    
    # if dat == 1:
    data4bc = data4bSig
    w3to4c = w3to4Sig
    predArrName_test = f'predSig_test_soft_{binInd}_nBins_{nBins}'
    predArrName_train = f'predSig_train_soft_{binInd}_nBins_{nBins}'
    data4b1, data4b2, w4b1, w4b2 = get_rand_pair(data4bc, binInd, weight = None)
    data3b1, data3b2, w3b1, w3b2 = get_rand_pair(data3b, binInd, weight = w3to4c)
    lab1, lab2 = np.concatenate((np.ones(len(w4b1)), np.zeros(len(w3b1)))) , np.concatenate((np.ones(len(w4b2)),np.zeros(len(w3b2))))
    data1, w1 = np.concatenate((data4b1, data3b1)), np.concatenate((w4b1, w3b1))
    data2, w2 = np.concatenate((data4b2, data3b2)), np.concatenate((w4b2, w3b2))
    # data2, w2 = np.concatenate((data4b2, data3b2)), np.concatenate((w4b2, w3b2))

    rand1 = np.random.choice(len(w1), size=len(w1), replace = False)
    rand2 = np.random.choice(len(w2), size=len(w2), replace = False)

    data_train, labels_train, weights_train = data1[rand1], lab1[rand1], w1[rand1]
    data_test , labels_test , weights_test  = data2[rand2], lab2[rand2], w2[rand2]

    # print(data_train,) 
    # print(labels_train,)
    # print(weights_train )                                                                         

    if not do_train:
        pred0 = (np.load(f'{predArrName_train}.npy')[~labels_train.astype(bool)]).reshape(-2) ## prob of bkg being data
        pred1 = (np.load(f'{predArrName_train}.npy')[labels_train.astype(bool)]).reshape(-2)  ## prob of dat being data
        wp0 = weights_train[~labels_train.astype(bool)]
        wp1 = weights_train[labels_train.astype(bool)]
        # b_pb_train = (2*pred0-1)#/(1-pred0)
        # d_pb_train = (2*pred1-1)#/(1-pred1)

        b_pb_train = (pred0)/(1-pred0)
        d_pb_train = (pred1)/(1-pred1)

        # b_pb_train = (pred0)/(1-pred0)
        # d_pb_train = (pred1)/(1-pred1)

        pred_test_b = (np.load(f'{predArrName_test}.npy')[~labels_test.astype(bool)]).reshape(-2) ## prob of bkg being data
        pred_test_d = (np.load(f'{predArrName_test}.npy')[labels_test.astype(bool)]).reshape(-2)  ## prob of dat being data
        wb_pb_test = weights_test[~labels_test.astype(bool)]
        wd_pb_test = weights_test[labels_test.astype(bool)]
        # b_pb_test = (2*pred_test_b-1)#/(1-pred_test_b)
        # d_pb_test = (2*pred_test_d-1)#/(1-pred_test_d)

        b_pb_test = (pred_test_b)/(1-pred_test_b)
        d_pb_test = (pred_test_d)/(1-pred_test_d)

        # b_pb_test = (pred_test_b)/(1-pred_test_b)
        # d_pb_test = (pred_test_d)/(1-pred_test_d)

        return b_pb_train, d_pb_train, wp0, wp1, b_pb_test, d_pb_test, wb_pb_test, wd_pb_test

    if do_train:
        model = keras.Sequential([
            layers.Input(shape=np.shape(data_train[0])),  
            # layers.Dense(64, activation=tf.keras.layers.LeakyReLU(alpha=0.01)), 
            # layers.Dropout(0.2),
            layers.Dense(32, activation=tf.keras.layers.LeakyReLU(alpha=0.01)),
            layers.Dropout(0.2),
            layers.Dense(16, activation=tf.keras.layers.LeakyReLU(alpha=0.01)),
            layers.Dropout(0.2),
            layers.Dense(4, activation=tf.keras.layers.LeakyReLU(alpha=0.01)),
            layers.Dense(1, activation='sigmoid')
        ])

        model.compile(optimizer='adam',
                    loss='binary_crossentropy',  
                    metrics=['accuracy'])

        model.summary()

        model.fit(data_train, labels_train, epochs=10, batch_size=200, 
                # validation_data=(data_valid, labels_valid, weights_valid),
                sample_weight = weights_train)

        pred_train = model.predict(data_train)
        pred_test = model.predict(data_test)

        del model
        K.clear_session()
        np.save(predArrName_train, pred_train)
        np.save(predArrName_test, pred_test)
        print('Done!')


# for binInd in range(10):
#     classify(binInd, do_train = True)
# exit()


ptest = PdfPages(f'test_score_prob_nBins_{nBins}.pdf')
ptrain = PdfPages(f'train_score_prob_nBins_{nBins}.pdf')

for binInd in range(10):
    xlab = 'p_sig/p_bkg'
    b_pb_train, d_pb_train, wp0, wp1, b_pb_test, d_pb_test, wb_pb_test, wd_pb_test = classify(binInd, do_train = False)
    ptrain = make_class_plot(ptrain, b_pb_train, d_pb_train, wp0, wp1, binInd = binInd, title = 'train' , xlab = xlab)
    ptest = make_class_plot(ptest, b_pb_test, d_pb_test, wb_pb_test, wd_pb_test, binInd = binInd, title = 'test', xlab = xlab)
ptest.close()
ptrain.close()    

exit()

def get_LRT_stat(pd_pb, count3b, count4b, norm ):
    return np.log(count3b/count4b) + (1/norm)*np.sum(np.log(pd_pb))

def get_LRT_stat_binned(pd_pb_cont, pd_pb_bounds, count3b, count4b, norm ):
    pd_pb_center = 0.5*(pd_pb_bounds[0:-1]+pd_pb_bounds[1:])
    return np.log(count3b/count4b) + (1/norm)*np.sum(pd_pb_cont*np.log(pd_pb_center))

def normal(x, mu = 0, std = 1):
    return np.exp(-0.5*((x - mu)/std)**2)/(std*np.sqrt(2*np.pi))


pprob = PdfPages(f'stats_sig_binned_nBins_{nBins}.pdf')
for binInd in range(10):
    b_pb_train, d_pb_train, wp0, wp1, b_pb_test, d_pb_test, wb_pb_test, wd_pb_test = classify(binInd, do_train = False)

    lo = np.min([np.min(b_pb_test), np.min(d_pb_test)])
    hi = np.max([np.max(b_pb_test), np.max(d_pb_test)])


    a = np.histogram(b_pb_test, range = (lo, hi), bins = 50, weights = wb_pb_test)
    b = np.histogram(d_pb_test, range = (lo, hi), bins = 50, weights = wd_pb_test);


    # stat_dat = get_LRT_stat(pd_pb=d_pb_test, count3b = np.sum(wb_pb_test), count4b = np.sum(wd_pb_test), norm = np.sum(wd_pb_test))
    # stat_bkg = get_LRT_stat(pd_pb=b_pb_test, count3b = np.sum(wb_pb_test), count4b = np.sum(wd_pb_test), norm = np.sum(wb_pb_test))

    stat_dat = get_LRT_stat_binned(a[0], a[1], count3b = np.sum(wb_pb_test), count4b = np.sum(wd_pb_test), norm = np.sum(wd_pb_test))
    stat_bkg = get_LRT_stat_binned(b[0], b[1], count3b = np.sum(wb_pb_test), count4b = np.sum(wd_pb_test), norm = np.sum(wb_pb_test))
    figf = plt.figure(figsize = (8,5))
    x = np.linspace(-5,5, 50)
    plt.plot(x, normal(x))
    plt.xlabel('LRT_stat_value')
    plt.title(f'binInd {binInd}')
    plt.axvline(stat_dat, ls = '--', color = 'g', label = f'stat_dat: {np.round(stat_dat,3)}')
    plt.axvline(stat_bkg, ls = '--', color = 'r', label = f'stat_bkg: {np.round(stat_bkg,3)}')
    plt.legend()
    pprob.savefig(figf)
    plt.close()

pprob.close()
    
    

