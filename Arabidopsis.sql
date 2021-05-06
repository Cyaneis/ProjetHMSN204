CREATE TABLE Genotype (
nom_genotype  varchar (30),
CONSTRAINT b PRIMARY KEY (Nom_genotype)
); 

\copy Genotype FROM 'genotype.csv' WITH CSV HEADER DELIMITER ',' QUOTE '"';


CREATE TYPE dm AS ENUM('Carence_Azote','Carence_Azote_Carbone','Standard','Saccharose', 'Microorganisme_PGPR') ; 

CREATE TABLE Milieu (
designation_milieu dm , 
nom_milieu  varchar (20), 
CONSTRAINT g PRIMARY KEY (nom_milieu)
); 

\copy Milieu FROM 'milieux.csv' WITH CSV HEADER DELIMITER ',' QUOTE '"';

CREATE TABLE Boite (
ID_boite int,
num_boite varchar (20), 
nom_milieu varchar (80),
CONSTRAINT f PRIMARY KEY (ID_boite), 
CONSTRAINT un FOREIGN KEY (nom_milieu) REFERENCES Milieu(nom_milieu)
); 

\copy Boite FROM 'boite.csv' WITH CSV HEADER DELIMITER ',' QUOTE '"';


CREATE TABLE Groupe (
ID_groupe int,
nb_plante_groupe int,
Groupe_Grouping  int, 
ID_boite  int, 
CONSTRAINT e PRIMARY KEY (ID_groupe),
CONSTRAINT deux FOREIGN KEY (ID_boite) REFERENCES Boite(ID_boite)
); 

\copy Groupe FROM 'groupe.csv' WITH CSV HEADER DELIMITER ',' QUOTE '"';

CREATE TABLE Caractere (
Nom varchar (50), 

CONSTRAINT c PRIMARY KEY (Nom)
);  

\copy Caractere FROM 'caracteres.csv' WITH CSV HEADER DELIMITER ',' QUOTE '"';

CREATE TABLE Plante (
ID_Plante varchar (10), 
ID_groupe int,
nom_genotype varchar (30),
CONSTRAINT a PRIMARY KEY (ID_Plante),
CONSTRAINT cinq FOREIGN KEY (ID_groupe) REFERENCES Groupe(ID_groupe)
); 

\copy Plante FROM 'plante.csv' WITH CSV HEADER DELIMITER ',' QUOTE '"';

CREATE TYPE uni AS ENUM('au','cm','cm_carre','none'); 
CREATE TABLE Mesure (
ID_mesure int,
Valeur real, 
unite uni, 
Nom_caractere varchar(50), 
ID_Plante varchar (10),
CONSTRAINT d PRIMARY KEY (ID_mesure),
CONSTRAINT trois FOREIGN KEY (Nom_caractere) REFERENCES Caractere(Nom),
CONSTRAINT quatre FOREIGN KEY (ID_plante) REFERENCES Plante(ID_plante)
); 



\copy Mesure FROM 'mesures.csv' WITH CSV HEADER DELIMITER ',' QUOTE '"';







