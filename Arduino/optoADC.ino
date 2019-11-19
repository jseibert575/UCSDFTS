/*
  Analog input, analog output, serial output
  Reads an analog input pin, maps the result to a range from 0 to 255 and uses
  the result to set the pulse width modulation (PWM) of an output pin which is optional.
  The program will print the results to the Serial Monitor or using the python script, can print to a CSV file to plot the data.
  This program was modified from existing tutorials from the Arduino website:
  
  https://www.arduino.cc/
  https://www.arduino.cc/en/Tutorial/AnalogInOutSerial
  https://www.arduino.cc/en/Tutorial/Calibration

*/

// These constants won't change. They're used to give names to the pins used:

const int analogInPin1 = A0;  // Analog input pin that the opto-interrupter is attached to
const int analogOutPin = 9;   // [[[optional]]] Digital output pin, If you want to output a digital signal using this pin.
const int analogInPin2 = A5;  // Analog input pin that reads the thermocouple value
int sensorValue = 0;          // voltage signal from the opto-interrupter 
                              // Arduino references from 5V as a default
                              // A value of 1024 corresponds to 5V signal


int outputValue = 0;          // value output to the PWM (analog out)
int sensorMin = 1023;
int sensorMax = 0;
unsigned long time;

int thermoCpl;                //thermocouple signal can be output if needed.

void setup(){
/*
  //OPTIONAL calibrate opto-interrupter for first 5 seconds
  //Measures max & min voltages read by opto-interrupter and calibrates readout values.
 
  while(millis() < 5000){
    sensorValue = analogRead(analogInPin1);
    if(sensorValue > sensorMax){
      sensorMax = sensorValue;
    }
    if(sensorValue < sensorMin){
      sensorMin = sensorValue;
    }
  }
*/

  // initialize serial communications at 115200 bps:
  //change this to change the output rate
    Serial.begin(115200);
}
void loop() {
  // read the analog in values:

  thermoCpl = analogRead(analogInPin2);
  sensorValue = analogRead(analogInPin1);
  
  // ADC steps - map it to the range of the analog out:
  outputValue = map(sensorValue, sensorMin, sensorMax, 0, 255);
  // change the analog out value:
  analogWrite(analogOutPin, outputValue);

  // print the data to the Serial Output:
  time = millis();                      // timestamp will be in milliseconds
  Serial.print(sensorValue);            // will print the analog signal being read
  Serial.print(" ");
  Serial.print(outputValue);            // will print the digital signal after ADC
  //Serial.print(" ");                  // Thermocouple data output
  //Serial.print(thermoCpl);
  Serial.print(" ");
  Serial.println(time);                 // prints time stamp and moves to the next line using println
  

 


  // Can add a delay to allow program to settle between reads just uncomment the delay
  // wait 1 milliseconds before the next loop for the analog-to-digital
  // converter to settle after the last reading:
  //delay(1);
}
