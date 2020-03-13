def objective(trial):
    clear_session()
    model = Sequential()
    model.add(Conv1D(filters=trial.suggest_int(
                         'filters1', 1, 1024),
                     kernel_size=trial.suggest_int(
                         'kernel_size1', 1, 42),
                     strides=trial.suggest_categorical(
                         'strides1', [1, 3, 5]),
                     activation=trial.suggest_categorical(
                         'activation1', ['linear', 'elu', 'relu', 'tanh']),
                     padding=trial.suggest_categorical(
                         'padding1', ['valid', 'same'])))
    model.add(MaxPooling1D(pool_size=trial.suggest_int(
                               'pool_size', 1, 42),
                           strides=trial.suggest_categorical(
                               'pool_strides', [1, 3, 5]),
                           padding=trial.suggest_categorical(
                               'pool_padding', ['valid', 'same'])))
    model.add(Conv1D(filters=trial.suggest_int(
                         'filters2', 1, 1024),
                     kernel_size=trial.suggest_int(
                         'kernel_size2', 1, 42),
                     strides=trial.suggest_categorical(
                         'strides2', [1, 3, 5]),
                     activation=trial.suggest_categorical(
                         'activation2', ['linear', 'elu', 'relu', 'tanh']),
                     padding=trial.suggest_categorical(
                         'padding2', ['valid', 'same'])))
    model.add(GlobalMaxPooling1D())
    model.add(Dropout(trial.suggest_discrete_uniform(
        'dropout', 0.0, 0.5, 0.25)))
    model.add(Dense(1, activation='sigmoid'))
    lr = trial.suggest_loguniform('lr', 1e-4, 1e-1)
    model.compile(loss='binary_crossentropy',
                  optimizer=Adam(learning_rate=lr),
                  metrics=['accuracy'])
    model.fit(X_train,
              y_train,
              validation_data=(X_test, y_test),
              shuffle=True,
              epochs=trial.suggest_categorical(
                  'epochs', [2, 4, 8, 16, 32, 64]),
              verbose=False)
    pred = model.predict_proba(X_test)
    score = roc_auc_score(y_test, pred)
    return score

