# ITU-R P.840-9 — Cloud & Fog Attenuation Model (Java)

This repository contains a **Java implementation of ITU-R P.840-9**, the International Telecommunication Union (ITU) Radiocommunication Sector recommendation for modeling attenuation due to clouds and fog in Earth–space communication systems.  

Official ITU sources:  
- [Official Recommendation ITU-R P.840-9 (08/2023)](https://www.itu.int/rec/R-REC-P.840-9-202308-I/en)  
- [ITU-R P.840-9 Dataset Readme File](https://www.itu.int/rec/R-REC-P.840Part15-0-202308-I/en)  
- [ITU P.840 Website](https://www.itu.int/rec/R-REC-P.840)  

It includes:  
- **Implementation** of the ITU-R P.840-9 prediction methods  
- **Testing** with unit tests  
- **Documentation** of the model and usage instructions  
- **Data files** (TXT maps from ITU recommendation and CSV test datasets)  

---

## 📂 Project Structure

```
ITU-R-P840-9/
├── data/                               # ITU dataset files (TXT format)
│   ├── annual/                         # Annual digital maps
│   ├── logNormalAnnual/                # Annual log-normal approximation digital maps
│   ├── month01/ ... month12/           # Monthly digital maps
│
├── docs/                               # Documentation (The generated Javadoc)
│
├── src/
│   ├── main/java/itu840/Itu840.java    # Main Java implementation
│   ├── test/java/itu840/...            # Unit tests
│   └── test/resources/                 # CSV test datasets for validation
│
├── pom.xml                             # Maven build file
```

---

## ⚙️ Features

The implementation covers all prediction methods of ITU-R P.840-9, including instantaneous, statistical, and log-normal approximation models.  
These features and more are explained in detail in the **generated Javadoc**.  

---

## 📊 Data Files

- **TXT files** (under `/data/`) → annual, monthly, and log-normal digital maps of climatological parameters.  
- **CSV files** (under `/src/test/resources/`) → used for automated testing and validation.  

The datasets and their usage are explained in detail in the official ITU documentation, its readme, and the implementation’s Javadoc. 

---

## 🚀 Usage

### Build
This is a Maven project. To compile:

```bash
mvn clean install
```

### Test
To execute the test suite:

```bash
mvn test
```

---

## 📜 License

This repository is released under the **MIT License**.  
You are free to use, modify, and distribute the code, provided that attribution is given.  

---

## ⚠️ Waiver

This implementation is based on the ITU-R P.840-9 recommendation and its associated datasets.  
The official documentation, datasets, and scientific definitions belong to the **International Telecommunication Union (ITU)**.  
This repository only provides an open-source **Java implementation** for educational and research purposes, and does not replace the official ITU sources.  

---

## 📖 Citation

If you use this implementation in research or projects, please cite as:

```
@software{sever_itu_r_p840_2025,
  author  = {Mehmet Sait Sever},
  title   = {Java Implementation of ITU-R P.840-9 (Cloud & Fog Attenuation Model)},
  year    = {2025},
  url     = {https://github.com/mss1091/ITU-R-P840-9}
}
```

---

## 📬 Contact

If you have any questions, suggestions, or recommendations, please do not hesitate to reach out to me at **severm21@itu.edu.tr**.  
