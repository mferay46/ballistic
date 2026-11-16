# -*- coding: utf-8 -*-
"""

Création de la pop up pour le projet ballistique
@author: FERAY.M
"""
import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QPushButton, QVBoxLayout,
    QHBoxLayout, QLabel
)
from PyQt5.QtGui import QPixmap

from PyQt5.QtCore import Qt


class TaskSelector(QWidget):

    def __init__(self):
        super().__init__()

        self.setWindowTitle("Choix du mode - Projet Balistique")
        self.selected_task = None      
        self.selected_mode = None       

        main_layout = QVBoxLayout()
        main_layout.setSpacing(20)
        self.setLayout(main_layout)


        title = QLabel("Sélection du calcul à effectuer")
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("font-size: 18px; font-weight: 600; margin-bottom:10px;")
        main_layout.addWidget(title)

        label_task = QLabel("1. Choisir la tâche :")
        main_layout.addWidget(label_task)

        task_layout = QHBoxLayout()
        main_layout.addLayout(task_layout)

        self.btn0 = QPushButton("0 : Modèle analytique")
        self.btn1 = QPushButton("1 : Numérique (drag/no_drag)")
        self.btn2 = QPushButton("2 : Trajectoires (cd, cl, T0)")
        self.btn3 = QPushButton("3 : Isocontours")

        for b in (self.btn0, self.btn1, self.btn2, self.btn3):
            b.setCheckable(True)
            task_layout.addWidget(b) 

        self.btn0.clicked.connect(lambda: self.set_task(0))
        self.btn1.clicked.connect(lambda: self.set_task(1))
        self.btn2.clicked.connect(lambda: self.set_task(2))
        self.btn3.clicked.connect(lambda: self.set_task(3))


        label_mode = QLabel("2.Uniquement pour le modèle numérique : Choisir le modèle drag/no drag et le mode de résolution :")
        main_layout.addWidget(label_mode)

        mode_layout = QHBoxLayout()
        main_layout.addLayout(mode_layout)

        self.btn_no = QPushButton("No_Drag réso : Odeint")
        self.btn_drag = QPushButton("Drag réso: Odeint")
        self.btn_no_ivp = QPushButton("No_Drag réso : IVP")
        self.btn_drag_ivp = QPushButton("Drag réso: IVP")

        for b in (self.btn_no, self.btn_drag,self.btn_no_ivp, self.btn_drag_ivp):
            b.setCheckable(True)
            mode_layout.addWidget(b)

        self.btn_no.clicked.connect(lambda: self.set_mode("no_drag"))
        self.btn_drag.clicked.connect(lambda: self.set_mode("drag"))
        self.btn_no_ivp.clicked.connect(lambda: self.set_mode("no_drag_ivp"))
        self.btn_drag_ivp.clicked.connect(lambda: self.set_mode("drag_ivp"))

        btn_ok = QPushButton("VALIDER")
        btn_ok.setFixedHeight(40)
        btn_ok.clicked.connect(self.close)
        main_layout.addWidget(btn_ok)

        btn_ok = QPushButton("Quitter le menu")
        btn_ok.setFixedHeight(30)
        btn_ok.clicked.connect(self.close)
        main_layout.addWidget(btn_ok)

  
    def set_task(self, t):
        """Active le bouton sélectionné et désactive les autres."""
        self.selected_task = t

        correspondance = {
            0: self.btn0,
            1: self.btn1,
            2: self.btn2,
            3: self.btn3,
        }

        for val, btn in correspondance.items():
            btn.setChecked(val == t)

        print(f"TASK sélectionné = {t}")

    def set_mode(self, mode):
        """Active le mode drag/no_drag."""
        self.selected_mode = mode

        self.btn_no.setChecked(mode == "no_drag")
        self.btn_drag.setChecked(mode == "drag")
        self.btn_drag_ivp.setChecked(mode == 'drag_ivp')
        self.btn_no_ivp.setChecked(mode == 'no_drag_ivp')

        print(f"MODE sélectionné = {mode}")




def start_popup():
    """
    Ouvre la fenêtre et retourne :
        - task (int)
        - mode (str : "drag" ou "no_drag")
    """
    app = QApplication(sys.argv)

    app.setStyleSheet("""
        QWidget {
            background-color: #46484d;
            color: #ffffff;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto;
            font-size: 13px;
        }

        QPushButton {
            background-color: #27272a;
            border: 1px solid #3f3f46;
            border-radius: 6px;
            padding: 8px 12px;
        }
        QPushButton:hover {
            background-color: #3f3f46;
        }
        QPushButton:checked {
            background-color: #4f46e5;
            border-color: #a5b4fc;
            color: white;
        }
    """)

    window = TaskSelector()
    window.resize(650, 300)
    window.show()

    app.exec_()

    return window.selected_task, window.selected_mode