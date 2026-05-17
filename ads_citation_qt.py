#!/usr/bin/env python3
"""Qt GUI wrapper for ADS citation metrics."""

from __future__ import annotations

import json
import sys
from typing import Any, Dict

from ads_citation_metrics import AdsClient, compute_metrics, extract_library_id

try:
    from PySide6.QtCore import QThread, Signal
    from PySide6.QtWidgets import (
        QApplication,
        QGridLayout,
        QGroupBox,
        QHBoxLayout,
        QLabel,
        QLineEdit,
        QMainWindow,
        QMessageBox,
        QPushButton,
        QScrollArea,
        QSpinBox,
        QTextEdit,
        QVBoxLayout,
        QWidget,
    )
except ImportError as exc:
    raise SystemExit(
        "PySide6 is required. Install it with: pip install PySide6"
    ) from exc


class Worker(QThread):
    success = Signal(dict)
    failed = Signal(str)

    def __init__(self, library_input: str, token: str, year_from: int, year_to: int) -> None:
        super().__init__()
        self.library_input = library_input
        self.token = token
        self.year_from = year_from
        self.year_to = year_to

    def run(self) -> None:
        try:
            library_id = extract_library_id(self.library_input)
            client = AdsClient(token=self.token)
            data = compute_metrics(
                client,
                library_id=library_id,
                year_from=self.year_from,
                year_to=self.year_to,
            )
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(str(exc))
            return
        self.success.emit(data)


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("ADS Citation Metrics")
        self.resize(900, 820)
        self._worker: Worker | None = None
        self._build_ui()

    def _build_ui(self) -> None:
        scroll_area = QScrollArea(self)
        scroll_area.setWidgetResizable(True)

        content = QWidget()
        root_layout = QVBoxLayout(content)
        root_layout.setContentsMargins(12, 12, 12, 12)
        root_layout.setSpacing(10)

        form_group = QGroupBox("Input")
        form_group.setObjectName("inputGroup")
        form_layout = QVBoxLayout(form_group)

        self.library_input = QLineEdit(
            "https://ui.adsabs.harvard.edu/public-libraries/7D7ggTUeS7aiE7dQzOHUkA"
        )
        self.library_input.setMinimumHeight(32)
        self.library_input.setMinimumWidth(450)
        self.token_input = QLineEdit()
        self.token_input.setEchoMode(QLineEdit.Password)
        self.token_input.setMinimumHeight(32)
        self.token_input.setMinimumWidth(450)

        library_label = QLabel("Library URL / ID")
        token_label = QLabel("ADS API Token")
        year_label = QLabel("Year Range")

        year_row = QWidget()
        year_layout = QHBoxLayout(year_row)
        year_layout.setContentsMargins(0, 0, 0, 0)
        year_layout.setSpacing(8)
        self.year_from_input = QSpinBox()
        self.year_from_input.setRange(1800, 2100)
        self.year_from_input.setValue(2021)
        self.year_to_input = QSpinBox()
        self.year_to_input.setRange(1800, 2100)
        self.year_to_input.setValue(2026)
        year_layout.addWidget(QLabel("From"))
        year_layout.addWidget(self.year_from_input)
        year_layout.addWidget(QLabel("To"))
        year_layout.addWidget(self.year_to_input)
        year_layout.addStretch(1)

        form_layout.addWidget(library_label)
        form_layout.addWidget(self.library_input)
        form_layout.addWidget(token_label)
        form_layout.addWidget(self.token_input)
        form_layout.addWidget(year_label)
        form_layout.addWidget(year_row)

        action_row = QWidget()
        action_layout = QHBoxLayout(action_row)
        action_layout.setContentsMargins(0, 0, 0, 0)
        action_layout.setSpacing(8)
        self.run_button = QPushButton("Compute")
        self.run_button.clicked.connect(self.on_run_clicked)
        self.status_label = QLabel("Ready.")
        action_layout.addWidget(self.run_button)
        action_layout.addWidget(self.status_label, 1)

        summary_group = QGroupBox("Summary")
        summary_group.setObjectName("summaryGroup")
        summary_layout = QGridLayout(summary_group)
        self.summary_labels: Dict[str, QLabel] = {}
        summary_fields = [
            ("papers_in_library", "Papers in library"),
            ("papers_in_year_range", "Papers in year range"),
            ("total_citations", "Total citations"),
            ("self_citations_detected", "Self citations detected"),
            ("total_citations_excluding_self", "Total excluding self"),
        ]
        for idx, (key, text) in enumerate(summary_fields):
            name = QLabel(text)
            value = QLabel("-")
            value.setStyleSheet("font-weight: 600;")
            summary_layout.addWidget(name, idx, 0)
            summary_layout.addWidget(value, idx, 1)
            self.summary_labels[key] = value

        top_group = QGroupBox("Top Paper (by citations excluding self)")
        top_group.setObjectName("topGroup")
        top_layout = QVBoxLayout(top_group)
        top_layout.setSpacing(6)
        self.top_title = QLineEdit()
        self.top_title.setReadOnly(True)
        self.top_title.setMinimumHeight(32)
        self.top_title.setMinimumWidth(450)
        self.top_bibcode = QLineEdit()
        self.top_bibcode.setReadOnly(True)
        self.top_bibcode.setMinimumHeight(32)
        self.top_bibcode.setMinimumWidth(450)
        self.top_year = QLineEdit()
        self.top_year.setReadOnly(True)
        self.top_year.setMinimumHeight(32)
        self.top_counts = QLineEdit()
        self.top_counts.setReadOnly(True)
        self.top_counts.setMinimumHeight(32)
        self.top_counts.setMinimumWidth(450)
        top_layout.addWidget(QLabel("Title"))
        top_layout.addWidget(self.top_title)
        top_layout.addWidget(QLabel("Bibcode"))
        top_layout.addWidget(self.top_bibcode)
        top_layout.addWidget(QLabel("Year"))
        top_layout.addWidget(self.top_year)
        top_layout.addWidget(QLabel("Counts"))
        top_layout.addWidget(self.top_counts)

        raw_group = QGroupBox("Raw JSON")
        raw_group.setObjectName("rawGroup")
        raw_layout = QVBoxLayout(raw_group)
        self.raw_output = QTextEdit()
        self.raw_output.setReadOnly(True)
        raw_layout.addWidget(self.raw_output)

        root_layout.addWidget(form_group)
        root_layout.addWidget(action_row)
        root_layout.addWidget(summary_group)
        root_layout.addWidget(top_group)
        root_layout.addWidget(raw_group, 1)
        scroll_area.setWidget(content)
        self.setCentralWidget(scroll_area)
        self.setStyleSheet(
            """
            QMainWindow {
                background: #f4f7fb;
            }
            QLabel {
                color: #1f2a44;
            }
            QGroupBox {
                font-weight: 600;
                border: 1px solid #c9d4ea;
                border-radius: 10px;
                margin-top: 12px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
                color: #16213d;
            }
            QGroupBox#inputGroup {
                background: #eaf3ff;
                border-color: #a9c7f5;
            }
            QGroupBox#summaryGroup {
                background: #eafaf1;
                border-color: #a6dfbe;
            }
            QGroupBox#topGroup {
                background: #fff8e8;
                border-color: #f0d38b;
            }
            QGroupBox#rawGroup {
                background: #f2efff;
                border-color: #c7bcf0;
            }
            QLineEdit, QSpinBox, QTextEdit {
                background: #ffffff;
                border: 1px solid #b8c7e6;
                border-radius: 6px;
                padding: 4px 8px;
                color: #1a2642;
            }
            QLineEdit:focus, QSpinBox:focus, QTextEdit:focus {
                border: 1px solid #4e89e8;
            }
            QPushButton {
                background: #3d7be0;
                color: white;
                border: 0;
                border-radius: 6px;
                padding: 8px 14px;
                font-weight: 600;
            }
            QPushButton:hover {
                background: #2f6ecf;
            }
            QPushButton:disabled {
                background: #9bb7e7;
            }
            """
        )

    def on_run_clicked(self) -> None:
        token = self.token_input.text().strip()
        library = self.library_input.text().strip()
        year_from = self.year_from_input.value()
        year_to = self.year_to_input.value()

        if not token:
            QMessageBox.warning(self, "Missing token", "Please input ADS API token.")
            return
        if not library:
            QMessageBox.warning(self, "Missing library", "Please input library URL/ID.")
            return
        if year_from > year_to:
            QMessageBox.warning(self, "Invalid years", "Year From must be <= Year To.")
            return

        self.run_button.setEnabled(False)
        self.status_label.setText("Computing... (this may take a while)")
        self._worker = Worker(library, token, year_from, year_to)
        self._worker.success.connect(self.on_success)
        self._worker.failed.connect(self.on_failed)
        self._worker.start()

    def on_success(self, data: Dict[str, Any]) -> None:
        self.run_button.setEnabled(True)
        self.status_label.setText("Done.")

        for key, label in self.summary_labels.items():
            value = data.get(key, "-")
            label.setText(f"{value:,}" if isinstance(value, int) else str(value))

        top = data.get("top_paper_by_citations_excluding_self") or {}
        self.top_title.setText(str(top.get("title", "-")))
        self.top_bibcode.setText(str(top.get("bibcode", "-")))
        self.top_year.setText(str(top.get("publication_year", "-")))

        total = top.get("total_citations", "-")
        self_c = top.get("self_citations", "-")
        non_self = top.get("citations_excluding_self", "-")
        self.top_counts.setText(
            f"Total: {total} | Self: {self_c} | Excluding self: {non_self}"
        )

        self.raw_output.setPlainText(json.dumps(data, indent=2, sort_keys=True))
        self._worker = None

    def on_failed(self, message: str) -> None:
        self.run_button.setEnabled(True)
        self.status_label.setText("Failed.")
        QMessageBox.critical(self, "Computation failed", message)
        self._worker = None


def main() -> int:
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
