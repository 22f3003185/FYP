# 📘 Guide to Using this Repository

---

## 📁 Folders in the Repo

| No. | Directory | Contents / Purpose |
|-----|-----------|--------------------|
| 1 | **Working** | This directory shall contain the finalized and validated codebase, comprising *solely* those features that have successfully completed all prescribed testing procedures. *No* new feature shall be introduced unless it has first undergone independent verification and subsequent integrated testing in conjunction with the materials contained in the **Experimentation Files** directory. |
| 2 | **Experimentation Files** | This directory shall serve as the controlled environment for the evaluation of candidate features that have already passed independent testing. *Only* those features that further demonstrate reliable and repeatable performance through extensive integrated testing herein shall be deemed *eligible for incorporation* into the finalized codebase within the **Working** directory. |
| 3 | **ROCOF Implementation** | This directory shall contain the implementation artifacts relating to ROCOF studies. The subdirectory **Generator Tripping – Modular** shall house the modular components designed to simulate the post-tripping behavior of a generator in a three-bus network. The standalone file **main.py** shall provide a consolidated implementation of the same functionality, together with the additional capability of plotting ROCOF and system frequency curves. |

---

## 🧪 Feature Testing and Promotion Procedure

### **1. Independent Verification Requirement**
Every proposed feature, modification, or enhancement shall undergo rigorous stand-alone testing to establish correctness, stability, and non-regression. Test evidence and outcomes shall be recorded.

### **2. Eligibility for Integrated Evaluation**
Following successful independent verification, the feature may proceed to integrated testing within the **Experimentation Files** directory of the same branch. This stage authorizes validation against existing and candidate components.

### **3. Standards During Integrated Testing**
Integrated evaluation shall confirm interoperability, repeatability, and acceptable system impact. Deficiencies shall mandate remediation and, where required, renewed independent qualification.

### **4. Approval for Promotion**
Only features meeting the defined criteria in both stages shall be eligible for transfer into the finalized implementations within the **Working** directory.

---

## 📄 Standalone Files

| No. | File | Contents / Purpose |
|-----|------|--------------------|
| 1 | **Bus Data.md** | This file shall serve as the authoritative reference for bus and load data associated with the IEEE 9-Bus system, presented in a concise and tabulated form for consistent reuse across the repository. |
| 2 | **FYP Dump.pdf** | This document shall preserve the collected background material, literature insights, and preliminary investigations assembled prior to the formal commencement of this project. |

---
