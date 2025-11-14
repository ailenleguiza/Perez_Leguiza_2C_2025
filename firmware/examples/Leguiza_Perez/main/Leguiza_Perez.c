/*! @mainpage Leguiza_Perez_EMG
 *
 * @section genDesc General Description
 *
 * Este proyecto ejemplifica el uso del módulo de comunicación 
 * Bluetooth Low Energy (BLE), junto con el cálculo de la FFT 
 * de una señal EMG adquirida desde un canal analógico.
 * Permite graficar en una aplicación móvil la FFT de la señal,
 * y detectar la fatiga muscular mediante el análisis espectral.
 *
 * @section changelog Changelog
 *
 * |   Date	    | Description                                    |
 * |:----------:|:-----------------------------------------------|
 * | 22/10/2025 | Código adaptado a EMG real con buffer circular |
 * | 28/10/2025 | Se agrega detección de fatiga mediante FFT     |
 * | 5/11/2025  | Se modifica codigo de la interfase             |
 *
 * @authors 
 * Florencia Ailen Leguiza Scandizzo  
 * María de los Ángeles Perez
 *
 */

/==================[inclusions]=============================================/
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

#include "led.h"
#include "ble_mcu.h"
#include "delay_mcu.h"

#include "fft.h"
#include "iir_filter.h"

#include "analog_io_mcu.h"
#include "timer_mcu.h"

/==================[macros and definitions]=================================/
#define CONFIG_BLINK_PERIOD 500
#define LED_BT	            LED_1
#define BUFFER_SIZE         512     // Ventana FFT
#define EMG_BUFFER_LEN      512    // Buffer circular
#define SAMPLE_FREQ	        512     // Hz
#define ADC_CH_EMG          CH1

// Parámetros de detección de fatiga
#define FATIGUE_THRESHOLD     0.15f   // 15% de descenso
#define BASELINE_WINDOWS      5       // Ventanas para calcular f_ref
#define CONSECUTIVE_WINDOWS   3       // Ventanas consecutivas necesarias

/==================[internal data definition]===============================/
static float emg_buffer[EMG_BUFFER_LEN];   // buffer circular para EMG
static uint16_t write_index = 0;
static uint16_t sample_count = 0;

static float emg_window[BUFFER_SIZE];       // ventana para FFT
static float emg_filt[BUFFER_SIZE];
static float emg_fft[BUFFER_SIZE/2];
static float emg_filt_fft[BUFFER_SIZE/2];
static float f[BUFFER_SIZE/2];

TaskHandle_t emg_task_handle = NULL;

/==================[variables para análisis de fatiga]======================/
static float f_ref = 0.0f;           // Frecuencia de referencia
static float f_mean = 0.0f;          // Frecuencia media actual
static float f_median = 0.0f;        // Frecuencia mediana actual
static float rms_value = 0.0f;       // RMS actual de la ventana
static int window_counter = 0;       // Contador de ventanas procesadas
static int fatigue_consecutive = 0;  // Ventanas consecutivas en fatiga
static bool fatigue_detected = false;

/**
 * @brief Función de callback que se ejecuta cuando se reciben datos por BLE.
 *
 * Si se recibe el carácter 'R', notifica la tarea EMGTask para procesar
 * una nueva ventana de datos EMG.  
 * Si se recibe el carácter 'B', se borra la gráfica de la FFT y
 * se limpian los vectores de datos asociados.
 *
 * @param data   Puntero al arreglo de bytes recibidos.
 * @param length Longitud del mensaje recibido.
 */
void read_data(uint8_t * data, uint8_t length){
    if(data[0] == 'R'){
        // Inicia una nueva medición
        xTaskNotifyGive(emg_task_handle);
    }
    else if(data[0] == 'B'){
        // Borra los gráficos en la app (usando el comando "C" que entiende Bluetooth Electronics)
        BleSendString("*HC*");   // limpia gráfico H (FFT cruda)
        vTaskDelay(pdMS_TO_TICKS(50)); // espera 50 ms para estabilidad
        // Envía mensaje de confirmación en el panel de texto
        BleSendString("*TLimpieza completada*\n");

    }
}


/**
 * @brief Escribe una muestra nueva en el buffer circular EMG.
 *
 * Guarda la muestra en la posición actual y avanza el índice de escritura.
 * Si el buffer aún no está lleno, incrementa el contador de muestras válidas.
 *
 * @param sample Valor leído del ADC (señal EMG).
 */
static void CircularBufferWrite(float sample){
    emg_buffer[write_index] = sample;
    write_index = (write_index + 1) % EMG_BUFFER_LEN;
    if(sample_count < EMG_BUFFER_LEN){
        sample_count++;
    }
}

/**
 * @brief Copia las últimas BUFFER_SIZE muestras desde el buffer circular
 *        hacia una ventana temporal para el procesamiento FFT.
 *
 * @param window Puntero al arreglo destino donde se copiarán las muestras.
 */
static void CircularBufferReadWindow(float *window){
    if(sample_count < BUFFER_SIZE) {
        return;
    }
    // Calcular inicio de las últimas BUFFER_SIZE muestras
    uint16_t start = (write_index - BUFFER_SIZE + EMG_BUFFER_LEN) % EMG_BUFFER_LEN;
    for(int i = 0; i < BUFFER_SIZE; i++){
        window[i] = emg_buffer[(start + i) % EMG_BUFFER_LEN];
    }
}

/==================[FUNCIONES AUXILIARES DE ANÁLISIS]======================/

/**
 * @brief Calcula el valor RMS (Root Mean Square) de un vector.
 *
 * Este valor representa la magnitud promedio de la señal,
 * y sirve como indicador de la activación muscular.
 *
 * @param data Puntero al vector de muestras.
 * @param len  Cantidad de elementos del vector.
 * @return Valor RMS.
 */
static float CalcRMS(float *data, int len){
    float sum = 0.0f;
    for(int i=0; i<len; i++) sum += data[i]*data[i];
    return sqrtf(sum / len);
}

/**
 * @brief Calcula la frecuencia media ponderada del espectro.
 *
 * Representa el centro de masa espectral y tiende a disminuir con la fatiga.
 *
 * @param fft_mag Vector de magnitudes espectrales.
 * @param freqs   Vector de frecuencias correspondientes.
 * @param len     Longitud de los vectores (BUFFER_SIZE/2).
 * @return Frecuencia media (Hz).
 */
static float CalcMeanFreq(float *fft_mag, float *freqs, int len){
    float num = 0.0f, den = 0.0f;
    for(int i=0; i<len; i++){
        num += freqs[i] * fft_mag[i];
        den += fft_mag[i];
    }
    return (den > 0) ? num / den : 0.0f;
}

/**
 * @brief Calcula la frecuencia mediana del espectro.
 *
 * Divide el área espectral en dos mitades iguales de energía.
 * Es una métrica robusta para estimar la fatiga muscular.
 *
 * @param fft_mag Vector de magnitudes espectrales.
 * @param freqs   Vector de frecuencias correspondientes.
 * @param len     Longitud de los vectores (BUFFER_SIZE/2).
 * @return Frecuencia mediana (Hz).
 */
static float CalcMedianFreq(float *fft_mag, float *freqs, int len){
    float total = 0.0f;
    for(int i=0; i<len; i++) total += fft_mag[i];
    float half = total / 2.0f;
    float accum = 0.0f;
    for(int i=0; i<len; i++){
        accum += fft_mag[i];
        if(accum >= half) return freqs[i];
    }
    return freqs[len - 1];
}

/**
 * @brief Tarea principal EMG: adquiere muestras, calcula FFT,
 *        obtiene métricas espectrales y detecta fatiga muscular.
 *
 * Se ejecuta cada vez que se recibe una notificación (por 'R' vía BLE).
 * Realiza el filtrado, cálculo de FFT, frecuencias características y
 * compara la frecuencia mediana con una referencia inicial.
 * Envía los resultados por Bluetooth Low Energy.
 *
 * @param pvParameter Parámetro de tarea (no utilizado).
 */
static void EMGTask(void *pvParameter){
    char msg_ble[128];
    static bool f_ref_established = false;  // ← Nueva bandera
    static float f_ref_accum = 0.0f;        // ← Acumulador temporal

    while(true){
        ulTaskNotifyTake(pdTRUE, portMAX_DELAY);

        // Tomar BUFFER_SIZE muestras desde buffer circular
        CircularBufferReadWindow(emg_window);

        // Filtros
        HiPassFilter(emg_window, emg_filt, BUFFER_SIZE);
        LowPassFilter(emg_filt, emg_filt, BUFFER_SIZE);

        // FFT
        FFTMagnitude(emg_window, emg_fft, BUFFER_SIZE);
        FFTMagnitude(emg_filt, emg_filt_fft, BUFFER_SIZE);
        FFTFrequency(SAMPLE_FREQ, BUFFER_SIZE, f);

        /==================[ANÁLISIS DE FATIGA Y ENVÍO BLE]=====================/
        window_counter++;

        /*==================[ENVÍO BLE DE LA FFT]=====================*/
        for(int i = 0; i < BUFFER_SIZE/2; i++){
            // Enviar datos FFT filtrada y cruda
            sprintf(msg_ble, "*HX%2.2fY%2.2f,X%2.2fY%2.2f*", 
                    f[i], emg_fft[i], f[i], emg_filt_fft[i]);
            BleSendString(msg_ble);
            vTaskDelay(pdMS_TO_TICKS(5));  // Evita saturar el stack BLE
        }


        // Calcular métricas de la ventana
        rms_value = CalcRMS(emg_filt, BUFFER_SIZE);
        f_mean    = CalcMeanFreq(emg_filt_fft, f, BUFFER_SIZE/2);
        f_median  = CalcMedianFreq(emg_filt_fft, f, BUFFER_SIZE/2);

        // Enviar métricas actuales por BLE

        sprintf(msg_ble, "*Tfmean: %.2fHz, fmed:%.2fHz, RMS:%.2f\n", f_mean, f_median, rms_value);
        BleSendString(msg_ble);
        // Enviar el número de ventana actual
        sprintf(msg_ble, "TVentana numero %d\n", window_counter);
        BleSendString(msg_ble);
        // --- FIN DE TU PETICIÓN ---


        // Calcular frecuencia de referencia (f_ref)
        if(window_counter <= BASELINE_WINDOWS){
            f_ref += f_median;
            if(window_counter == BASELINE_WINDOWS){
                f_ref /= BASELINE_WINDOWS;
                sprintf(msg_ble, "*Tf_ref establecida: %.2f Hz*\n", f_ref);
                printf("Frecuencia referencia: %.2f Hz\n", f_ref);
                BleSendString(msg_ble);
            }
        }
        // Luego del baseline, comparar con umbral
        else if(f_ref > 0.0f){
            float drop = (f_ref - f_median) / f_ref;
            sprintf(msg_ble, "*T Drop: %.2f\n", drop);
            BleSendString(msg_ble);

            // Evaluar fatiga
            if (drop > FATIGUE_THRESHOLD) {
                fatigue_consecutive++;

                if (!fatigue_detected) {
                    sprintf(msg_ble, "*TCaida detectada (%.1f%% > %.1f%%)*\n", drop * 100.0f, FATIGUE_THRESHOLD * 100.0f);
                    BleSendString(msg_ble);
                }

                if (fatigue_consecutive >= CONSECUTIVE_WINDOWS && !fatigue_detected) {
                    fatigue_detected = true;
                    sprintf(msg_ble, "*TFATIGA DETECTADA: (decae%.1f%% respecto a referencia)*\n", drop * 100.0f);
                    BleSendString(msg_ble);
                }

            }        
            BleSendString(msg_ble);

        /==================[FASE DE BASELINE]=====================/
        if(!f_ref_established){
            if(window_counter <= BASELINE_WINDOWS){
                f_ref_accum += f_median;
                sprintf(msg_ble, "TBase %.0d/%.0d: f_med=%.2f\n", window_counter, BASELINE_WINDOWS, f_median);
                BleSendString(msg_ble);
            }

            if(window_counter == BASELINE_WINDOWS){
                f_ref = f_ref_accum / BASELINE_WINDOWS;
                f_ref_established = true;
                sprintf(msg_ble, "Tf_ref establecida: %.2f Hz\n", f_ref);
                BleSendString(msg_ble);
                printf("Frecuencia de referencia establecida: %.2f Hz\n", f_ref);
            }
        }

        /==================[FASE DE DETECCIÓN DE FATIGA]=====================/
        else {
            if(f_ref > 0.0f){
                float drop = (f_ref - f_median) / f_ref;
                sprintf(msg_ble, "*T Drop: %.2f\n", drop);
                BleSendString(msg_ble);

                // Evaluar fatiga
                if (drop > FATIGUE_THRESHOLD) {
                    fatigue_consecutive++;

                    if (!fatigue_detected) {
                        sprintf(msg_ble, "TCaida detectada (%.1f%% > %.1f%%)\n", drop * 100.0f, FATIGUE_THRESHOLD * 100.0f);
                        BleSendString(msg_ble);
                    }

                    if (fatigue_consecutive >= CONSECUTIVE_WINDOWS && !fatigue_detected) {
                        fatigue_detected = true;
                        sprintf(msg_ble, "TFATIGA DETECTADA (decae %.1f%% respecto a ref)\n", drop * 100.0f);
                        BleSendString(msg_ble);
                    }
                } else {
                    fatigue_consecutive = 0; // se resetea si se recupera
                }
            }

        }
    }
}


/**
 * @brief Rutina de interrupción del temporizador.
 *
 * Se ejecuta a la frecuencia de muestreo definida y adquiere
 * una muestra analógica del canal EMG, almacenándola en el buffer circular.
 *
 * @param param Parámetro del temporizador (no utilizado).
 */
void EMG_TimerISR(void *param){
    uint16_t adc_val;
    AnalogInputReadSingle(ADC_CH_EMG, &adc_val);
    CircularBufferWrite((float)adc_val);
}

/**
 * @brief Función principal de la aplicación.
 *
 * Inicializa los periféricos: LEDs, BLE, ADC, temporizador y FFT.
 * Crea la tarea EMG encargada del procesamiento de señal y detección de fatiga.
 * El bucle principal actualiza el estado del LED según el estado del BLE.
 */
void app_main(void){
    // Inicializaciones
    LedsInit();
    FFTInit();
    LowPassInit(SAMPLE_FREQ, 30, ORDER_2);
    HiPassInit(SAMPLE_FREQ, 1, ORDER_2);

    // BLE
    ble_config_t ble_configuration = {
        "ESP_EMG",
        read_data
    };
    BleInit(&ble_configuration);

    // ADC
    analog_input_config_t adc_config = {
        .input = ADC_CH_EMG,
        .mode = ADC_SINGLE,
        .func_p = NULL,
        .param_p = NULL,
        .sample_frec = SAMPLE_FREQ
    };
    AnalogInputInit(&adc_config);

    // Timer para muestreo EMG
    timer_config_t emg_timer = {
        .timer = TIMER_A,
        .period = 1000000 / SAMPLE_FREQ, // en us
        .func_p = EMG_TimerISR,
        .param_p = NULL
    };
    TimerInit(&emg_timer);
    TimerStart(TIMER_A);

    // Tarea EMG
    xTaskCreate(&EMGTask, "EMG", 4096, NULL, 5, &emg_task_handle);

    // Loop principal: indica estado BLE
    while(1){
        vTaskDelay(CONFIG_BLINK_PERIOD / portTICK_PERIOD_MS);
        switch(BleStatus()){
            case BLE_OFF:
                LedOff(LED_BT);
            break;
            case BLE_DISCONNECTED:
                LedToggle(LED_BT);
            break;
            case BLE_CONNECTED:
                LedOn(LED_BT);
            break;
        }
    }
}