import subprocess
import time
from datetime import datetime
import numpy as np
import pyvisa

class AdbBridge:
    def __init__(self, ip_puerto):
        self.ip_puerto = ip_puerto

    def conectar(self):
        """Conecta ADB al dispositivo."""
        print(f"[ADB] Conectando a {self.ip_puerto}...")
        resultado = subprocess.run(["adb", "connect", self.ip_puerto], capture_output=True, text=True)
        print(resultado.stdout.strip())

    def tomar_captura(self, nombre_archivo=None):
        """Toma una captura de pantalla del dispositivo."""
        if not nombre_archivo:
            nombre_archivo = datetime.now().strftime("screenshot_%Y%m%d_%H%M%S.png")

        print(f"[ADB] Tomando captura y guardando como {nombre_archivo}...")
        with open(nombre_archivo, "wb") as f:
            subprocess.run(["adb", "-s", self.ip_puerto, "exec-out", "screencap", "-p"], stdout=f)

    def esta_conectado(self):
        """Verifica si el dispositivo sigue conectado."""
        resultado = subprocess.run(["adb", "devices"], capture_output=True, text=True)
        return self.ip_puerto.split(":")[0] in resultado.stdout

    def touch_center(self):
        subprocess.run(["adb", "shell", "input", "tap", str(512), str(300)])
        
    def adb_disconnect(self):
        subprocess.run(["adb", "disconnect", self.ip_puerto], capture_output=True, text=True)

class DHO804:
    def __init__(self,ip,adbport,scpi_port=None):
        self.ip = ip
        # scpi_port (p. ej. 5555): SCPI por socket TCP en vez de VXI-11. Sirve cuando
        # el servicio VXI-11 del osciloscopio queda colgado (no responde ni *IDN?)
        # tras cortar una transferencia a mitad de camino; el socket sigue andando.
        self.resource = (f'TCPIP::{ip}::{scpi_port}::SOCKET' if scpi_port
                         else f'TCPIP::{ip}::INSTR')
        self.adb_bridge = AdbBridge(ip + ':' + adbport)
        self.rm = pyvisa.ResourceManager()    

    def connect(self):
        """Conecta interfaz ADB y SCPI.

        Reintenta el open VXI-11: el WiFi del osciloscopio tiene power-save y la
        radio puede estar dormida en el momento de conectar. El timeout VISA
        generoso tolera los despertares de radio durante el barrido.
        """
        self.adb_bridge.conectar()
        last = None
        for attempt in range(4):
            try:
                self.scpi_bridge = self.rm.open_resource(
                    self.resource, open_timeout=8000)
                self.scpi_bridge.write_termination = '\n'
                self.scpi_bridge.read_termination = '\n'
                self.scpi_bridge.timeout = 8000
                return
            except Exception as e:
                last = e
                time.sleep(1.5)
        raise last

    def take_screenshot(self):
        if (self.adb_bridge.esta_conectado()):
            self.adb_bridge.tomar_captura()
        else:
            return 0

    def auto_set(self):
        self.scpi_bridge.write(":AUToset")
        self.adb_bridge.touch_center()

    def measure_avg(self, ch: int):
        if self._validate_channel(ch):
            return float(self.scpi_bridge.query(f':MEASure:VAVG? CHANnel{ch}'))
        else:
            pass

    def measure_amp(self, ch: int):
        if self._validate_channel(ch):
            return float(self.scpi_bridge.query(f':MEASure:VAMP? CHANnel{ch}'))
        else:
            pass
    def measure_rms(self, ch: int):
        if self._validate_channel(ch):
            return float(self.scpi_bridge.query(f':MEASure:VRMS? CHANnel{ch}'))
        else:
            pass

    def measure_vpp(self, ch: int):
        if self._validate_channel(ch):
            return float(self.scpi_bridge.query(f':MEASure:VPP? CHANnel{ch}'))

    def measure_phase(self, ch_a: int, ch_b: int):
        """Fase de ch_a respecto de ch_b en grados (flanco de subida vs subida)."""
        if self._validate_channel(ch_a) and self._validate_channel(ch_b):
            return float(self.scpi_bridge.query(
                f':MEASure:ITEM? RRPHase,CHANnel{ch_a},CHANnel{ch_b}'))

    def set_channel_display(self, ch: int, on: bool):
        if self._validate_channel(ch):
            self.scpi_bridge.write(f':CHANnel{ch}:DISPlay {"ON" if on else "OFF"}')

    def set_channel_scale(self, ch: int, scale: float):
        """Escala vertical en V/div."""
        if self._validate_channel(ch):
            self.scpi_bridge.write(f':CHANnel{ch}:SCALe {scale}')

    def get_channel_scale(self, ch: int):
        if self._validate_channel(ch):
            return float(self.scpi_bridge.query(f':CHANnel{ch}:SCALe?'))

    def set_channel_offset(self, ch: int, offset: float):
        """Offset vertical en V."""
        if self._validate_channel(ch):
            self.scpi_bridge.write(f':CHANnel{ch}:OFFSet {offset}')

    def set_channel_coupling(self, ch: int, coupling: str):
        """Acople de entrada del canal: AC, DC o GND."""
        if self._validate_channel(ch):
            self.scpi_bridge.write(f':CHANnel{ch}:COUPling {coupling}')

    def set_channel_probe(self, ch: int, ratio: float):
        """Relación de atenuación de la sonda (p.ej. 1 o 10).

        El DHO800 valida PROBe contra una lista discreta con match estricto:
        acepta "1"/"10" pero rechaza "1.0"/"10.0" (-222 Data out of range). Se
        formatea con :g para colapsar el .0 de los enteros sin romper los ratios
        fraccionarios (0.1, 0.001).
        """
        if self._validate_channel(ch):
            self.scpi_bridge.write(f':CHANnel{ch}:PROBe {ratio:g}')

    def set_timebase(self, scale: float):
        """Escala horizontal (base de tiempo principal) en s/div."""
        self.scpi_bridge.write(f':TIMebase:MAIN:SCALe {scale}')

    def set_acquire_average(self, count: int):
        """Promediado de adquisición (count = potencia de 2). count<=1 -> NORMal."""
        if count and count > 1:
            self.scpi_bridge.write(':ACQuire:TYPE AVERages')
            self.scpi_bridge.write(f':ACQuire:AVERages {count}')
        else:
            self.scpi_bridge.write(':ACQuire:TYPE NORMal')

    def set_trigger_edge(self, source: int, level: float = 0.0, slope: str = "POSitive"):
        """Trigger por flanco sobre un canal: referencia estable para promediar y fase."""
        if self._validate_channel(source):
            self.scpi_bridge.write(':TRIGger:MODE EDGE')
            self.scpi_bridge.write(':TRIGger:SWEep AUTO')
            self.scpi_bridge.write(f':TRIGger:EDGE:SOURce CHANnel{source}')
            self.scpi_bridge.write(f':TRIGger:EDGE:SLOPe {slope}')
            self.scpi_bridge.write(f':TRIGger:EDGE:LEVel {level}')

    def set_memory_depth(self, mdepth: str):
        """Profundidad de memoria: AUTO, 1k, 10k, 100k, 1M, 5M o 10M (con 2 canales, hasta 5M)."""
        self.scpi_bridge.write(f':ACQuire:MDEPth {mdepth}')

    def run(self):
        self.scpi_bridge.write(':RUN')

    def single(self, timeout: float = 30.0, force_after: float = None):
        """Dispara una adquisición única y espera a que termine (estado STOP).

        Primero espera a que el estado SALGA de STOP: justo después de :SINGle el
        osciloscopio todavía puede contestar el STOP de la adquisición anterior, y
        leer en ese momento devuelve la memoria vieja (o vacía). Con force_after (s),
        si el disparo no engancha en ese tiempo, lo fuerza con :TFORce.
        """
        b = self.scpi_bridge
        b.write(':SINGle')
        t0 = time.time()
        while b.query(':TRIGger:STATus?').strip() == 'STOP':
            if time.time() - t0 > 2.0:
                break                      # ya terminó (adquisición muy rápida)
            time.sleep(0.02)
        forced = False
        while b.query(':TRIGger:STATus?').strip() != 'STOP':
            elapsed = time.time() - t0
            if force_after is not None and not forced and elapsed > force_after:
                b.write(':TFORce')
                forced = True
            if elapsed > timeout:
                raise TimeoutError('El osciloscopio no completó la adquisición única')
            time.sleep(0.05)

    def read_waveform(self, ch: int, chunk: int = 100_000):
        """Lee la memoria interna del canal (modo RAW, formato WORD) en tensión.

        El osciloscopio tiene que estar en STOP (ver single()). Devuelve (v, dt):
        array de tensiones en V y el intervalo de muestreo en s. Se lee por
        bloques porque la memoria puede tener millones de puntos.
        Tensión = (dato - YORigin - YREFerence) * YINCrement (guía de programación,
        sección 3.28). La guía no dice el orden de bytes del WORD: se asume
        little-endian; thd.py lo verifica contra :MEASure:VPP?.
        """
        if not self._validate_channel(ch):
            raise ValueError(f'Canal inválido: {ch}')
        b = self.scpi_bridge
        b.write(f':WAVeform:SOURce CHANnel{ch}')
        b.write(':WAVeform:MODE RAW')
        b.write(':WAVeform:FORMat WORD')
        # Esperar a que aplique modo y formato, y leer la escala de una sola vez con
        # PREamble: consultando YINC/YOR/YREF por separado justo después de cambiar
        # de modo, el DHO804 a veces devuelve valores viejos (amplitud mal escalada).
        # La cantidad de puntos sale de la memoria: el campo "points" del PREamble es
        # el tamaño de la ventana STARt-STOP que quedó de la lectura anterior.
        npts = int(float(b.query(':ACQuire:MDEPth?')))
        b.write(':WAVeform:STARt 1')
        b.write(f':WAVeform:STOP {min(chunk, npts)}')
        b.query('*OPC?')
        pre = b.query(':WAVeform:PREamble?').split(',')
        # format, type, points, count, xinc, xorigin, xref, yinc, yorigin, yref
        dt = float(pre[4])
        yinc, yor, yref = float(pre[7]), float(pre[8]), float(pre[9])
        old_timeout = b.timeout
        b.timeout = 30000
        try:
            data = []
            for start in range(1, npts + 1, chunk):
                stop = min(start + chunk - 1, npts)
                b.write(f':WAVeform:STARt {start}')
                b.write(f':WAVeform:STOP {stop}')
                block = b.query_binary_values(
                    ':WAVeform:DATA?', datatype='H', is_big_endian=False,
                    container=np.array)
                # Con memoria menor que el bloque, el DHO804 no siempre respeta STOP
                # y devuelve el bloque entero: los puntos de más son basura.
                n_req = stop - start + 1
                if len(block) < n_req:
                    raise IOError(f'Bloque incompleto: {len(block)} de {n_req} puntos')
                data.append(block[:n_req])
        finally:
            b.timeout = old_timeout
        raw = np.concatenate(data).astype(np.float64)
        return (raw - yor - yref) * yinc, dt

    def _validate_channel(self,ch:int):
        return (ch > 0 and ch <=5)

    def close(self):
        self.scpi_bridge.close()
        self.rm.close()
        self.adb_bridge.adb_disconnect()
