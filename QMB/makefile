# Variables de configuración (CamelCase obligatorio)
OutputSource = QMB_Source.txt
FileExtensions = -name "*.m" -o -name "*.wl"

# Objetos virtuales
.PHONY: all clean $(OutputSource)

# Objetivo principal
all: $(OutputSource)

$(OutputSource):
	@echo "Iniciando la extracción del código fuente de QMB..."
	@find . -type f \( $(FileExtensions) \) | while read -r File; do \
		echo "=== FILENAME: QMB/$${File#./} ==="; \
		cat "$$File"; \
		echo "-e "; \
		echo ""; \
	done > $(OutputSource)
	@echo "¡Archivo $(OutputSource) generado exitosamente en el directorio!"

clean:
	@echo "Eliminando $(OutputSource)..."
	rm -f $(OutputSource)
