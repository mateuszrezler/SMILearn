from tensorflow.keras.callbacks import Callback


class EpochCounter(Callback):

    def on_epoch_begin(self, epoch, logs):
        print(f'\rrunning epoch #{epoch+1}', end='')

    def on_epoch_end(self, epoch, logs):
        print('\r', end='')

