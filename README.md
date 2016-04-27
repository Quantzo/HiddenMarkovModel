# HiddenMarkovModel
##Program functions
```
x <- train(Markov, observationTraining[1:15000], 100)
evaluateDataSet(testSet, x$hmm)
```
##Results
#####Uniform growth model
> - Uniform starting probability
- Test data set length : 2 000
- Train data set length : 3 130
- Train data sets observation number : 518 550
- Used observartions : 15 000

| Helix states number | Strand states number | Coil statas number | Result |
| ------------------- | -------------------- | ------------------ | ------ |
| 1 | 1 | 1 | 0.3401149 |
| 2 | 2 | 2 | 0.3805153 |
| 3 | 3 | 3 | 0.3383858 |
| 4 | 4 | 4 |  |
| 5 | 5 | 5 |  |
| 6 | 6 | 6 |  |
| 7 | 7 | 7 |  |
| 8 | 8 | 8 |  |
