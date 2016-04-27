import sqlite3
import json
from pathlib import Path

class Protein:
    def __init__(self, id, firstStructure, secondaryStructure):
        self.id = id
        self.firstStructure = firstStructure
        self.secondaryStructure = self.__alterSecondaryStructure__(secondaryStructure)

    def __replaceMultipleChars__(self, str, chars, newValue):
        for c in chars:
            if c in str:
                str = str.replace(c, newValue)
        return str

    def __alterSecondaryStructure__(self, str):
        helixes = "GI"
        strands = "E"
        coils = "-ST"
        str = self.__replaceMultipleChars__(str, helixes, 'H')
        str = self.__replaceMultipleChars__(str, strands, 'B')
        str = self.__replaceMultipleChars__(str, coils, 'C')
        return str


                    



def addToDatabase(object, cur):
    cur.execute("insert into proteins values(:id, :firstStructure, :secondaryStructure)", object)

def createDataBase():
    connection = sqlite3.connect('proteins.db')
    cursor = connection.cursor()
    cursor.execute('''CREATE TABLE proteins (id text, firstStructure text, secondaryStructure text)''')
    connection.commit()
    connection.close()

def step(line,file):
    return Protein(line.rstrip(), file.readline().rstrip(), file.readline().rstrip())

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

