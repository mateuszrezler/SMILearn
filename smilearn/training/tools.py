from tensorflow.keras.callbacks import Callback
from scipy.stats import ttest_ind


class EpochCounter(Callback):

    def on_epoch_begin(self, epoch, logs):
        print(f'\rrunning epoch #{epoch+1}', end='')

    def on_epoch_end(self, epoch, logs):
        print('\r', end='')


def compare_tstudent(list1, list2, p=0.05):
    stat, pt = ttest_ind(list1, list2)
    print('The means are', 'not'*int(p < pt), 'different.')

