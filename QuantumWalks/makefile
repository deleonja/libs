# Variables
OUTPUT = QuantumWalks_Source.txt
EXTENSIONS = -name "*.m" -o -name "*.wl"

# Forzamos a que estos objetivos siempre se ejecuten
.PHONY: all clean $(OUTPUT)

# El objetivo por defecto ahora es generar el archivo
all: $(OUTPUT)

$(OUTPUT):
	@echo "Recreando $(OUTPUT) desde cero..."
	@find . -type f \( $(EXTENSIONS) \) | while read -r file; do \
		echo "=== FILENAME: QuantumWalks/$${file#./} ==="; \
		cat "$$file"; \
		echo -e "\n"; \
	done > $(OUTPUT)
	@echo "¡Archivo $(OUTPUT) generado exitosamente!"

clean:
	rm -f $(OUTPUT)
