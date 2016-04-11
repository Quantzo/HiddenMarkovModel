import sqlite3
import json
from pathlib import Path

class Protein:
    def __init__(self, id, firstStructure, secondaryStructure):
        self.id = id
        self.firstStructure = firstStructure
        self.secondaryStructure = secondaryStructure

def addToDatabase(object, cur):
    cur.execute("insert into proteins values(:id, :firstStructure, :secondaryStructure)", object)

def createDataBase():
    connection = sqlite3.connect('proteins.db')
    cursor = connection.cursor()
    cursor.execute('''CREATE TABLE proteins (id text, firstStructure text, secondaryStructure text)''')
    connection.commit()
    connection.close()

def step(line,file):
    return Protein(line, file.readline(), file.readline())

def readFile(fileName,cursor, jsonList):
    with fileName.open() as file:
        x = True
        while(x):
            line = file.readline()
            if(line == ''):
                return
            else:
                protein = step(line, file)
                addToDatabase(protein.__dict__, cursor)
                jsonList.append(protein.__dict__)


def checkIfDataBaseExist():
    dataBase = Path('proteins.db')
    return dataBase.exists()

def readFiles():

    if(not checkIfDataBaseExist()):
        createDataBase()
    dir = Path('Proteins')
    files = list(dir.glob('**/*.txt'))
    
    connection = sqlite3.connect('proteins.db')
    cursor = connection.cursor()

    jsonList = list()

    for file in files:
        readFile(file,cursor, jsonList)

    connection.commit()
    connection.close()


    with open("jsonProteins.json", 'w') as jsonFile:
        json.dump(jsonList, jsonFile)





readFiles()

