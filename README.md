# HiddenMarkovModel
##Program functions
```
x <- train(Markov, observationTraining[1:15000], 100)
evaluateDataSet(testSet, x$hmm)
```
##Tests

> - Test 1 : Uniform Growth test, where n = <1,16>

> | Helix states number | Strand states number | Coil states number | 
| ------------------- | -------------------- | ------------------ |
| 13 | 13 | 13 |

> - Test 2 : Uniform Growth test with constant number of Helixes states, where n = <7,15>

> | Helix states number | Strand states number | Coil states number | 
| ------------------- | -------------------- | ------------------ |
| 13 | 14 | 14 |

> - Test 2 : Uniform Growth test with constant number of Helixes and Strands, where n = <10,18>

> | Helix states number | Strand states number | Coil states number | 
| ------------------- | -------------------- | ------------------ |
| 13 | 14 | n |

##Results
#####Uniform growth model
> - Uniform starting probability
- Test data set length : 2 000
- Train data set length : 3 130
- Train data sets observation number : 518 550
- Used observartions : 15 000

| Helix states number | Strand states number | Coil states number | Result |
| ------------------- | -------------------- | ------------------ | ------ |
| 1 | 1 | 1 | 0.3401149 |
| 2 | 2 | 2 | 0.3805153 |
| 3 | 3 | 3 | 0.3383858 |
| 4 | 4 | 4 | 0.35377308707124 |
| 5 | 5 | 5 | 0.225953748253919 |
| 6 | 6 | 6 | 0.433931398416887 |
| 7 | 7 | 7 | 0.340114853329194 |
